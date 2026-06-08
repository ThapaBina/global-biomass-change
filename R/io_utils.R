read_tile()
read_reference()
save_model()
load_model()
save_raster()
ensure_dir()
#---
# =========================
# IO UTILITIES
# =========================

read_raster <- function(path) {
  terra::rast(path)
}

write_raster <- function(x, filename, overwrite = TRUE) {
  terra::writeRaster(x, filename, overwrite = overwrite)
}

save_model <- function(model, path) {
  saveRDS(model, path)
}

load_model <- function(path) {
  readRDS(path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

get_tile_id <- function(file_path) {
  tools::file_path_sans_ext(basename(file_path))
}
