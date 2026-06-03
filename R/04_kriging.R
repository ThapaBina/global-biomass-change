# =========================
# 04_kriging.R
# =========================

library(terra)
library(gstat)

source("R/utils.R")

rf_resid_3km <- rast("outputs/rf_resid_3km.tif")

rf_resid_utm <- project_to_local_utm(rf_resid_3km)

df <- as.data.frame(rf_resid_utm, xy = TRUE, na.rm = TRUE)

vgm_emp <- variogram(dAGB ~ 1, ~x+y, data = df)

vgm_fit <- fit.variogram(
  vgm_emp,
  vgm(psill = 0.07, model = "Sph", range = 200000, nugget = 0.1)
)

gs <- gstat(NULL, "dAGB", dAGB ~ 1, df, locations = ~x+y, model = vgm_fit)

krig <- interpolate(rf_resid_utm, gs)

writeRaster(krig, "outputs/krig_3km.tif", overwrite = TRUE)