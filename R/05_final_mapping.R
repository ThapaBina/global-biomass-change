# =========================
# 05_final_mapping.R
# =========================

library(terra)

rf_pred_30m <- readRDS("outputs/rf_pred_30m.rds")
krig_3km    <- rast("outputs/krig_3km.tif")
ba_3km      <- readRDS("outputs/ba_3km.rds")

fact <- 100  # adjust if needed

krig_30m <- resample(krig_3km, rf_pred_30m)

final_30m <- rf_pred_30m + krig_30m

writeRaster(final_30m, "outputs/final_30m.tif", overwrite = TRUE)

final_3km <- aggregate(final_30m, fact, mean)
writeRaster(final_3km, "outputs/final_3km.tif", overwrite = TRUE)