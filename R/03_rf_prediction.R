# =========================
# 03_rf_prediction.R
# =========================

library(terra)
library(caret)

rf_model <- readRDS("outputs/rf_model.rds")
pred_3km <- readRDS("outputs/pred_3km.rds")
ba_3km   <- readRDS("outputs/ba_3km.rds")

# 3km prediction
rf_pred_3km <- predict(pred_3km, rf_model$finalModel, type = "quantiles", quantiles = 0.5)

rf_resid_3km <- ba_3km - rf_pred_3km

writeRaster(rf_pred_3km, "outputs/rf_pred_3km.tif", overwrite = TRUE)
writeRaster(rf_resid_3km, "outputs/rf_resid_3km.tif", overwrite = TRUE)