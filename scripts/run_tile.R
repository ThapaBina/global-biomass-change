# =========================================================
# RF-ATPK GLOBAL BIOMASS CHANGE WORKFLOW
# Single tile execution script
# =========================================================

suppressPackageStartupMessages({
  library(terra)
  library(gstat)
  library(caret)
  library(ranger)
  library(ggplot2)
  library(viridisLite)
  library(parallel)
})

# =========================
# LOAD UTILITIES
# =========================

source("R/io_utils.R")
source("R/raster_utils.R")
source("R/modeling_utils.R")
source("R/mapping_utils.R")

# =========================
# MAIN FUNCTION
# =========================

run_tile <- function(tile_id,
                     tile_file,
                     reference_file,
                     output_dir) {

  cat("\n=========================================\n")
  cat("Processing tile:", tile_id, "\n")
  cat("=========================================\n")

  ensure_dir(output_dir)
  
  # managing the temporary directory for rasters:
  terra_tmp_dir <- file.path(output_dir, "tmp_terra")
  
  dir.create(terra_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  terra_tmp_dir <- normalizePath(terra_tmp_dir, winslash = "/", mustWork = FALSE)
  
  terraOptions(tempdir = terra_tmp_dir, memfrac = 0.5,  progress = 1)
  
  Sys.setenv(TMPDIR = terra_tmp_dir, GDAL_CACHEMAX = "2048"  )
  

  # -------------------------
  # 1. READ INPUT DATA
  # -------------------------

  fine_raster <- read_raster(tile_file)         # 30 m predictor (GLAD)
  coarse_ref  <- read_raster(reference_file)    # 3 km AGB reference

  coarse_ref <- crop_to_reference(coarse_ref, fine_raster)
  names(coarse_ref) <- "dAGB"
  
  #
  fine_raster
  coarse_ref
  
  # crop
  fact <- compute_aggregation_factor(fine = fine_raster, coarse=coarse_ref)
  cat("Aggregation factor:", fact, "\n")
  cat("Random Forest Modeling:\n")
  
  # -------------------------
  # 2. BUILD PREDICTORS
  # -------------------------

  predictors_3km <- create_predictors(fine = fine_raster, coarse = coarse_ref, fact = fact)

  # -------------------------
  # 3. TRAIN RANDOM FOREST
  # -------------------------

  train_df <- build_training_df(response = coarse_ref,predictors = predictors_3km )

  rf_model <- train_rf(train_df)

  save_model(rf_model,file.path(output_dir, paste0(tile_id, "_rf_model.rds")))

  # -------------------------
  # 4. RF PREDICTION (3 km)
  # -------------------------
  cat("RF Prediction (3 km)\n")
  
  # rf_pred_3km <- terra::predict(predictors_3km, rf_model$finalModel,
  #                               fun = function(model, data) {predict(model, data)$predictions})
  rf_pred_3km<- terra::predict(predictors_3km, rf_model$finalModel,
                               type = "quantiles", quantiles = c(0.5))

  # -------------------------
  # 5. RESIDUALS (3 km)
  # -------------------------
  cat("RF residual computing ....\n")

  rf_residuals_3km <- compute_residuals(coarse_ref,  rf_pred_3km)
  
  rf_residuals_3km_file <- file.path(output_dir, paste0(tile_id, "_RF_residuals_3km.tif"))
  
  write_raster(rf_residuals_3km, rf_residuals_3km_file)

  # -------------------------
  # 6. VARIOGRAM FITTING
  # -------------------------
  cat("Variogram and Kriging Residual Modeling:\n")

  vg_obj <- fit_variogram(rf_residuals_3km)

  save_model(vg_obj$model, file.path(output_dir, paste0(tile_id, "_variogram.rds")) )

  # -------------------------
  # 7. KRIGING (RESIDUALS)
  # -------------------------

  krig_result <- krige_residuals(rf_residuals_3km, vg_obj )
  
  # VARIANCE / UNCERTAINTY MAP
  rf_krig_residuals <- krig_result[[2]]

  rf_krig_file <- file.path(output_dir,paste0(tile_id, "_krig_variance_3km.tif"))

  write_raster(rf_krig_residuals, rf_krig_file)
  
  # extract the Kriging residuals:
  Krig_pred_3km <- krig_result[[1]]
  
  ## remove files to save space:
  rm(rf_residuals_3km, krig_result, rf_krig_residuals)
  
  # -------------------------
  # 8. RF PREDICTION (30 m)
  # -------------------------
  cat("Random Forest Model Prediction: 30 m\n")
  # this part takes long time 
  rf_pred_30m = predict_rf_raster_fine_resolution(rf_model, fine_raster)
  
  # ------------------------------------------
  # 9. COMBINE RF + KRIG PREDICTION (30 m)
  # -----------------------------------------
  cat("Final Prediction: RF + Krig\n")
  # resample Kriging residual to 30 m 
  Krig_pred_30m <- resample(Krig_pred_3km, rf_pred_30m, method = "bilinear")
  
  # add (RF prediction + Krig Prediction)
  final_30m <- rf_pred_30m + Krig_pred_30m
  
  final_30m_file <- file.path(output_dir,paste0(tile_id, "_RFATPKpred_30m.tif"))

  write_raster(final_30m, final_30m_file, overwrite = TRUE,
               wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")))
  
  # ----------------------------------------
  # 10. UPSCALE COMBINED 3 km → 30 m
  # ----------------------------------------
  final_3km <- aggregate(final_30m, fact, mean, na.rm = TRUE)
  final_3km <- resample(final_3km, coarse_ref, method = "bilinear")
  
  final_3km_file <- file.path(output_dir,paste0(tile_id, "_RFATPKpred_3km.tif"))
  
  writeRaster(final_3km, final_3km_file, overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")))

  # -------------------------------------------------
  # 11. MASS PRESERVATION CORRECTION
  # -------------------------------------------------
  cat("Mass Preservation Correction\n")
  ## 30 m raster
  final_30m_cor = mass_preservation_correction(final_3km, coarse_ref, final_30m)
  
  final_30m_file_cor <- file.path(output_dir,paste0(tile_id, "_RFATPKpred_30m_cor.tif"))
  
  write_raster(final_30m_cor, final_30m_file_cor, overwrite = TRUE,
               wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")))
  
  
  # Aggregated corrected 3 km
  final_3km_cor <- aggregate(final_30m_cor, fact, mean, na.rm = TRUE)
  final_3km_cor <- resample(final_3km_cor, coarse_ref, method = "bilinear")
  
  final_3km_file_cor <- file.path(output_dir,paste0(tile_id, "_RFATPKpred_3km_cor.tif"))
  
  write_raster(final_3km_cor, final_3km_file_cor, overwrite = TRUE,
               wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")))
  
  
  # -------------------------
  # 12. SCATTER PLOT
  # -------------------------

  cat("Plotting....\n")
  
  df_plot <- as.data.frame(c(coarse_ref, final_3km, final_3km_cor),xy = TRUE, na.rm = TRUE)

  colnames(df_plot)[3:5] <- c("obs", "pred", "corr")

  p1 <- plot_scatter(df_plot, "pred", "obs",
                     paste0(tile_id, " RF + Kriging"))

  ggsave(filename = file.path(output_dir, paste0(tile_id, "_scatter.png")),
         plot = p1, width = 6,  height = 5 )

  # -------------------------
  # DONE
  # -------------------------

  cat("\nTile completed:", tile_id, "\n")

  return(list(
    rf_model = rf_model,
    variogram = vg_obj$model,
    final_map = final_30m_cor
  ))
}

# =========================
# CLI ENTRY POINT
# =========================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 4) {

  run_tile(
    tile_id = args[1],
    tile_file = args[2],
    reference_file = args[3],
    output_dir = args[4]
  )

} else {

  cat("Usage:\n")
  cat("Rscript run_tile.R <tile_id> <tile_file> <reference_file> <output_dir>\n")
}
