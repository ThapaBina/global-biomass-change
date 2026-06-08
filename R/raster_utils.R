# =========================
# RASTER UTILITIES
# =========================

crop_to_reference <- function(r, ref) {
  terra::crop(r, terra::ext(ref))
}

aggregate_raster <- function(r, fact, fun = mean) {
  terra::aggregate(r, fact = fact, fun = fun, na.rm = TRUE)
}

aggregate_quantiles <- function(r, fact, probs = c(0.1, 0.9)) {

  qfun <- function(x, ...) {
    quantile(x, probs = probs, na.rm = TRUE)
  }

  terra::aggregate(r, fact = fact, fun = qfun)
}

resample_to_reference <- function(r, ref, method = "bilinear") {
  terra::resample(r, ref, method = method)
}

project_to_utm <- function(r) {

  e <- terra::ext(r)

  lon <- (e[1] + e[2]) / 2
  lat <- (e[3] + e[4]) / 2

  zone <- floor((lon + 180) / 6) + 1

  epsg <- if (lat >= 0) 32600 + zone else 32700 + zone

  terra::project(r, paste0("EPSG:", epsg))
}

match_resolution <- function(source, target) {
  terra::resample(source, target, method = "bilinear")
}

