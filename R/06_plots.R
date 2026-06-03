# =========================
# 06_plots.R
# =========================

library(terra)
library(ggplot2)
library(viridisLite)

final_30m <- rast("outputs/final_30m.tif")

png("figures/final_map.png", width = 800, height = 600)
plot(final_30m, col = viridis(100))
dev.off()