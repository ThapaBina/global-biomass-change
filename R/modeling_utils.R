# =========================================================
# MODELING UTILITIES (RF + KRIGING + CORRECTION)
# =========================================================

# -------------------------
# RANDOM FOREST
# -------------------------

# =========================================================
# CREATE PREDICTORS
# =========================================================

# =========================================================
# CREATE PREDICTORS AT COARSE RESOLUTION
# =========================================================
#
# fine   = fine-resolution raster (e.g., GLAD 30 m)
# coarse = reference raster (e.g., dAGB 3 km)
# fact   = aggregation factor
#
# Returns:
#   SpatRaster with:
#     tch_mean
#     tch_sd
#     tch_q10
#     tch_q90
#
# =========================================================

create_predictors <- function(
    fine,
    coarse,
    fact
) {
  
  # ----------------------------------------
  # Quantile function
  # ----------------------------------------
  
  q_fun <- function(x, ...) {
    stats::quantile(
      x,
      probs = c(0.1, 0.9),
      na.rm = TRUE
    )
  }
  
  # ----------------------------------------
  # Aggregate fine-resolution raster
  # ----------------------------------------
  
  tch_mean <- terra::aggregate(
    fine,
    fact = fact,
    fun = mean,
    na.rm = TRUE
  )
  
  tch_sd <- terra::aggregate(
    fine,
    fact = fact,
    fun = sd,
    na.rm = TRUE
  )
  
  tch_q <- terra::aggregate(
    fine,
    fact = fact,
    fun = q_fun
  )
  
  # ----------------------------------------
  # Rename quantile layers
  # ----------------------------------------
  
  names(tch_q) <- c(
    "tch_q10",
    "tch_q90"
  )
  
  # ----------------------------------------
  # Match reference grid
  # ----------------------------------------
  
  tch_mean <- terra::resample(
    tch_mean,
    coarse,
    method = "bilinear"
  )
  
  tch_sd <- terra::resample(
    tch_sd,
    coarse,
    method = "bilinear"
  )
  
  tch_q <- terra::resample(
    tch_q,
    coarse,
    method = "near"
  )
  
  # ----------------------------------------
  # Combine predictors
  # ----------------------------------------
  
  predictors <- c(
    tch_mean,
    tch_sd,
    tch_q
  )
  
  names(predictors) <- c(
    "tch_mean",
    "tch_sd",
    "tch_q10",
    "tch_q90"
  )
  
  return(predictors)
}


build_training_df <- function(response, predictors) {
  as.data.frame(c(response, predictors), xy = TRUE, na.rm = TRUE)
}

train_rf <- function(df) {

  ctrl <- caret::trainControl(
    method = "cv",
    number = 5
  )

  grid <- expand.grid(
    mtry = 2,
    splitrule = "variance",
    min.node.size = 5
  )

  caret::train(
    dAGB ~ tch_mean + tch_sd + tch_q10 + tch_q90,
    data = df,
    method = "ranger",
    trControl = ctrl,
    tuneGrid = grid,
    num.trees = 1000,
    importance = "permutation"
  )
}

predict_rf_raster <- function(model, predictors) {
  terra::predict(
    predictors,
    model$finalModel,
    type = "response"
  )
}

# =========================================================
# CREATE 30 m PREDICTION STACK
# =========================================================

create_prediction_stack_30m <- function(fine_raster) {
  
  r_stack <- c(
    fine_raster,
    terra::rast(fine_raster),
    terra::rast(fine_raster),
    terra::rast(fine_raster)
  )
  
  names(r_stack) <- c(
    "tch_mean",
    "tch_sd",
    "tch_q10",
    "tch_q90"
  )
  
  return(r_stack)
}



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

  # vg_fit <- gstat::fit.variogram(
  #   vg_emp,
  #   gstat::vgm(
  #     model = "Sph",
  #     psill = 0.07,
  #     range = 200000,
  #     nugget = 0.1
  #   )
  # )
  vg_fit <- gstat::fit.variogram(vg_emp, model = gstat::vgm('Sph'))

  list(
    empirical = vg_emp,
    model = vg_fit,
    data = df
  )
}

# -------------------------
# KRIGING
# -------------------------
# vgm_fit_cp_fast <- gstat(NULL,"dAGB", dAGB ~ 1, resid_df, locations = ~x+y,
#                          model = vgm_fit,
#                          nmax = 500,
#                          maxdist = vgm_fit$range[2])

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



combine_rf_kriging <- function(rf_pred, krig_obj) {
  rf_pred + krig_obj[[1]]
}

# -------------------------
# MASS PRESERVATION
# -------------------------

compute_correction_factor <- function(final_3km, reference_3km) {
  reference_3km / final_3km
}

# apply_correction <- function(final_30m, correction_3km) {
# 
#   correction_30m <- terra::resample(
#     correction_3km,
#     final_30m,
#     method = "bilinear"
#   )
# 
#   final_30m * correction_30m
# }

apply_correction <- function(
    final_30m,
    correction_3km
) {
  
  fact <- round(
    terra::res(correction_3km)[1] /
      terra::res(final_30m)[1]
  )
  
  correction_30m <- terra::disagg(
    correction_3km,
    fact
  )
  
  correction_30m <- terra::resample(
    correction_30m,
    final_30m,
    method = "bilinear"
  )
  
  final_30m * correction_30m
}
