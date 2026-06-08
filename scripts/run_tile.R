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

  # -------------------------
  # 1. READ INPUT DATA
  # -------------------------

  fine_raster <- read_raster(tile_file)         # 30 m predictor (GLAD)
  coarse_ref  <- read_raster(reference_file)    # 3 km AGB reference

  coarse_ref <- crop_to_reference(coarse_ref, fine_raster)
  names(coarse_ref) <- "dAGB"

  # -------------------------
  # 2. BUILD PREDICTORS
  # -------------------------

  predictors_3km <- create_predictors(
    fine = fine_raster,
    coarse = coarse_ref
  )

  # -------------------------
  # 3. TRAIN RANDOM FOREST
  # -------------------------

  train_df <- build_training_df(
    response = coarse_ref,
    predictors = predictors_3km
  )

  rf_model <- train_rf(train_df)

  save_model(
    rf_model,
    file.path(output_dir, paste0(tile_id, "_rf_model.rds"))
  )

  cat("RF model trained\n")

  # -------------------------
  # 4. RF PREDICTION (3 km)
  # -------------------------

  rf_pred_3km <- predict_rf_raster(
    rf_model,
    predictors_3km
  )

  rf_pred_3km_file <- file.path(output_dir,
                                paste0(tile_id, "_RF_pred_3km.tif"))

  write_raster(rf_pred_3km, rf_pred_3km_file)

  # -------------------------
  # 5. RESIDUALS (3 km)
  # -------------------------

  rf_residuals_3km <- compute_residuals(
    coarse_ref,
    rf_pred_3km
  )

  # -------------------------
  # 6. VARIOGRAM FITTING
  # -------------------------

  vg_obj <- fit_variogram(rf_residuals_3km)

  save_model(
    vg_obj$model,
    file.path(output_dir, paste0(tile_id, "_variogram.rds"))
  )

  cat("Variogram fitted\n")

  # -------------------------
  # 7. KRIGING (RESIDUALS)
  # -------------------------

  krig_result <- krige_residuals(
    rf_residuals_3km,
    vg_obj
  )

  rf_krig_residuals <- krig_result[[1]]

  rf_krig_file <- file.path(output_dir,
                            paste0(tile_id, "_krig_residuals_3km.tif"))

  write_raster(rf_krig_residuals, rf_krig_file)

  # -------------------------
  # 8. COMBINE RF + KRIGING
  # -------------------------

  rf_plus_krig_3km <- combine_rf_kriging(
    rf_pred_3km,
    krig_result
  )

  rf_plus_krig_file <- file.path(output_dir,
                                 paste0(tile_id, "_RF_Krig_3km.tif"))

  write_raster(rf_plus_krig_3km, rf_plus_krig_file)

  # -------------------------
  # 9. RF PREDICTION (30 m)
  predictor_stack_30m <- create_prediction_stack_30m(
    fine_raster
  )
  
  rf_pred_30m <- predict_rf_raster(
    rf_model,
    predictor_stack_30m
  )
  
  rf_pred_30m_file <- file.path(output_dir,
                                paste0(tile_id, "_RF_pred_30m.tif"))

  write_raster(rf_pred_30m, rf_pred_30m_file)

  # -------------------------
  # 10. UPSCALE COMBINED 3 km → 30 m
  # -------------------------

  rf_krig_30m <- match_resolution(
    rf_plus_krig_3km,
    rf_pred_30m
  )

  # -------------------------
  # 11. FINAL PREDICTION (30 m)
  # -------------------------

  final_30m <- rf_pred_30m + rf_krig_30m

  final_30m_file <- file.path(output_dir,
                              paste0(tile_id, "_FINAL_30m.tif"))

  write_raster(final_30m, final_30m_file)

  # -------------------------
  # 12. AGGREGATE BACK TO 3 km
  # -------------------------

  fact <- round(res(coarse_ref)[1] / res(fine_raster)[1])

  final_3km <- aggregate_raster(final_30m, fact, mean)

  final_3km <- resample_to_reference(final_3km, coarse_ref)

  final_3km_file <- file.path(output_dir,
                              paste0(tile_id, "_FINAL_3km.tif"))

  write_raster(final_3km, final_3km_file)

  # -------------------------
  # 13. MASS PRESERVATION CORRECTION
  # -------------------------

  correction_factor <- coarse_ref / final_3km

  final_30m_corr <- apply_correction(
    final_30m,
    correction_factor
  )

  final_corr_file <- file.path(output_dir,
                               paste0(tile_id, "_FINAL_30m_CORR.tif"))

  write_raster(final_30m_corr, final_corr_file)

  # -------------------------
  # 14. VARIANCE / UNCERTAINTY MAP
  # -------------------------

  krig_variance <- sqrt(krig_result[[2]])

  var_file <- file.path(output_dir,
                        paste0(tile_id, "_KRIG_VARIANCE_3km.tif"))

  write_raster(krig_variance, var_file)

  # -------------------------
  # 15. SCATTER PLOT
  # -------------------------

  df_plot <- as.data.frame(
    c(coarse_ref, final_3km, final_3km),
    xy = TRUE,
    na.rm = TRUE
  )

  colnames(df_plot)[3:5] <- c("obs", "pred", "corr")

  p1 <- plot_scatter(df_plot, "pred", "obs",
                     paste0(tile_id, " RF + Kriging"))

  ggsave(
    filename = file.path(output_dir,
                         paste0(tile_id, "_scatter.png")),
    plot = p1,
    width = 6,
    height = 5
  )

  # -------------------------
  # DONE
  # -------------------------

  cat("\nTile completed:", tile_id, "\n")

  return(list(
    rf_model = rf_model,
    variogram = vg_obj$model,
    final_map = final_30m_corr
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
