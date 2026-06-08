# =========================
# utils.R
# =========================

library(terra)
library(sf)
library(gstat)
library(dplyr)

# -------------------------
# UTM projection helper
# -------------------------
project_to_local_utm <- function(rast) {
  
  ext <- ext(rast)
  lon_cent <- (ext[1] + ext[2]) / 2
  lat_cent <- (ext[3] + ext[4]) / 2
  
  zone <- floor((lon_cent + 180) / 6) + 1
  epsg <- if (lat_cent >= 0) 32600 + zone else 32700 + zone
  
  project(rast, paste0("EPSG:", epsg))
}

# -------------------------
# quantile function
# -------------------------
q_fun <- function(x) {
  quantile(x, probs = c(0.1, 0.9), na.rm = TRUE)
}