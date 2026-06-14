# =========================================================
# MODELING UTILITIES (RF + KRIGING + CORRECTION)
# =========================================================

# fine   = fine-resolution raster (e.g., GLAD 30 m)
# coarse = reference raster (e.g., dAGB 3 km)
# fact = aggregation factor

# =========================================================
# Create RF predictors at 3 km
# =========================================================
create_predictors <- function(fine, coarse, fact) {
  
  # ----------------------------------------
  # Quantile function
  # ----------------------------------------
  q_fun <- function(x, ...) {stats::quantile( x, probs = c(0.1, 0.9), na.rm = TRUE) }
  
  # ----------------------------------------
  # Aggregate fine-resolution raster
  # ----------------------------------------
  tch_mean <- terra::aggregate(fine,fact = fact,fun = mean,  na.rm = TRUE )
  
  tch_sd <- terra::aggregate(fine,fact = fact,fun = sd, na.rm = TRUE)
  
  tch_q <- terra::aggregate(fine,  fact = fact,  fun = q_fun)
  
  # ----------------------------------------
  # Match reference grid
  # ----------------------------------------
  tch_mean <- terra::resample(tch_mean,  coarse,  method = "bilinear" )
  
  tch_sd <- terra::resample(tch_sd,  coarse,  method = "bilinear" )
  
  tch_q <- terra::resample(tch_q,coarse,   method = "near"  )
  
  # ----------------------------------------
  # Combine predictors
  # ----------------------------------------
  predictors <- c(tch_mean, tch_sd,  tch_q)
  names(predictors) <- c( "tch_mean",  "tch_sd","tch_q10","tch_q90" )
  
  return(predictors)
}

# =========================================================
# Create dataframe for RF modeling
# =========================================================
build_training_df <- function(response, predictors) {
  as.data.frame(c(response, predictors), xy = TRUE, na.rm = TRUE)
}

train_rf <- function(df) {
  ctrl <- caret::trainControl(method = "cv", number = 5 )

  grid <- expand.grid(mtry = 2,  splitrule = "variance",  min.node.size = 5 )

  caret::train(dAGB ~ tch_mean + tch_sd + tch_q10 + tch_q90, data = df, 
               method = "ranger", trControl = ctrl, tuneGrid = grid,
               num.trees = 1000, importance = "permutation", quantreg=TRUE, num.threads=10 )
}


# =========================================================
# CREATE 30 m PREDICTION STACK
# =========================================================

create_prediction_stack_30m <- function(fine_raster) {
  
  r_stack <- c(fine_raster, terra::rast(fine_raster),  terra::rast(fine_raster),terra::rast(fine_raster))
  names(r_stack) <- c("tch_mean", "tch_sd",  "tch_q10",  "tch_q90")
  return(r_stack)
}

# =========================================================
# Compute RF Residuals for Spatial Analysis: Vario + Krig
# =========================================================
compute_residuals <- function(obs, pred) {
  obs - pred
}

# -------------------------
# VARIOGRAM
# -------------------------
fit_variogram <- function(residual_raster) {
  
  r_utm <- project_to_utm(residual_raster)
  
  df <- as.data.frame(r_utm, xy = TRUE, na.rm = TRUE)
  
  vg_emp <- gstat::variogram(dAGB ~ 1, ~x + y, data = df)
  
  vg_fit <- gstat::fit.variogram(vg_emp, model = gstat::vgm("Sph"))
  
  # plot empirical + fitted variogram
  print(plot(vg_emp, vg_fit))
  
  invisible(list(
    empirical = vg_emp,
    model = vg_fit,
    data = df
  ))
}

# -------------------------
# KRIGING
# -------------------------

krige_residuals <- function(residual_raster,  vg_obj) {
  
  r_utm <- project_to_utm(residual_raster)
  
  df <- vg_obj$data
  
  #sp::coordinates(df) <- ~x+y
  
  g <- gstat(NULL, 'dAGB',
    formula = dAGB ~ 1,
    data = df, locations = ~x+y,
    model = vg_obj$model,
    nmax = 500,
    maxdist = vg_obj$model$range[2]
  )
  
  krig_utm <- terra::interpolate(r_utm,  g)
  
  terra::project(krig_utm, residual_raster)
}


# -------------------------
# RF PREDICTION (30 m)
# -------------------------

predict_rf_raster_fine_resolution <- function(model, data, tile_id, output_dir) {
  
  #library(terra)
  # dir.create('./Output/terra_tmp')
  # terraOptions(tempdir = "./Output/terra_tmp")
  
  rf_model <- model$finalModel
  
  # ---------------------------
  # 1. Create tiles
  # ---------------------------
  grid_r <- rast(ext(data), nrows = 25, ncols = 25, crs = crs(data))
  tiles <- as.polygons(grid_r)
  
  cat("Number of tiles:", nrow(tiles), "\n")
  
  # ---------------------------
  # 2. Temp folder
  # ---------------------------
  tmp_dir <- file.path(output_dir,'tmp')
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------
  # 3. Prediction function (IMPORTANT)
  # ---------------------------
  pred_fun <- function(mod, df) {
    predict(mod, data = df)$predictions
  }
  
  # ---------------------------
  # 4. Loop over tiles (FIXED): seq_along(tiles)
  # ---------------------------
  for (i in seq_along(tiles)[c(15,51:55)]) {
    
    r_sub <- crop(data, tiles[i])
    r_stack <- c(r_sub, rast(r_sub), rast(r_sub),rast(r_sub))
    names(r_stack) = c("tch_mean", "tch_sd", "tch_q10", "tch_q90")
    
    terra::predict(
      r_stack,
      rf_model,
      fun = pred_fun,
      filename = file.path(tmp_dir, sprintf("pred_tile_%04d.tif", i)),
      overwrite = TRUE,
      wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=NONE"))
    )
    
    rm(r_sub, r_stack)
    gc()   # 🔥 critical fix for crash after few tiles
  }
  
  # ---------------------------
  # 5. Mosaic safely
  # ---------------------------
  pred_files <- list.files(tmp_dir, pattern = "\\.tif$", full.names = TRUE)
  
  mosaic_r <- rast(pred_files[1])
  
  for (f in pred_files[-1]) {
    mosaic_r <- mosaic(mosaic_r, rast(f))
  }
  
  # ---------------------------
  # 6. Write final output
  # ---------------------------
  #out_file <- "./Outputs/rf_prediction_final.tif"
  out_file <- file.path(output_dir, paste0(tile_id, "_RF_pred_30m.tif"))
  
  writeRaster(mosaic_r, out_file, overwrite = TRUE,
              wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW")) )
  
  unlink(tmp_dir, recursive = TRUE)
  
  cat("Prediction done!!!\n")
  
  return(mosaic_r)
}



# ------------------------------------------------------------------
# MASS PRESERVATION
# -------------------------------------------------------------------

mass_preservation_correction <- function(final_3km, coarse_ref, final_30m){
  correction_factor = coarse_ref / final_3km
  fact <- compute_aggregation_factor(fine = final_30m, coarse=coarse_ref)
  correction_factor_30m <- disagg(correction_factor, fact)
  correction_factor_30m <- resample(correction_factor_30m, final_30m, method = "bilinear")
  final_30m_corr <- final_30m * correction_factor_30m
  return(final_30m_corr)
}

