# =========================================================
# MODELING UTILITIES (RF + KRIGING + CORRECTION)
# =========================================================

# -------------------------
# RANDOM FOREST
# -------------------------

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

  vg_fit <- gstat::fit.variogram(
    vg_emp,
    gstat::vgm(
      model = "Sph",
      psill = 0.07,
      range = 200000,
      nugget = 0.1
    )
  )

  list(
    empirical = vg_emp,
    model = vg_fit,
    data = df
  )
}

# -------------------------
# KRIGING
# -------------------------

krige_residuals <- function(residual_raster, vg_obj) {

  df <- vg_obj$data
  vgm <- vg_obj$model

  g <- gstat::gstat(
    NULL,
    "resid",
    dAGB ~ 1,
    data = df,
    model = vgm
  )

  terra::interpolate(residual_raster, g)
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

apply_correction <- function(final_30m, correction_3km) {

  correction_30m <- terra::resample(
    correction_3km,
    final_30m,
    method = "bilinear"
  )

  final_30m * correction_30m
}
