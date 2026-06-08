# =========================================================
# MAPPING UTILITIES
# =========================================================

plot_raster <- function(r, title = "") {
  terra::plot(
    r,
    col = viridisLite::viridis(256),
    main = title
  )
}

plot_scatter <- function(df, x, y, title = "") {

  fit <- lm(df[[y]] ~ df[[x]])

  r2 <- summary(fit)$r.squared
  rmse <- sqrt(mean((df[[y]] - df[[x]])^2, na.rm = TRUE))

  ggplot2::ggplot(df, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", color = "red") +
    ggplot2::labs(title = title) +
    ggplot2::annotate(
      "text",
      x = Inf, y = -Inf,
      label = paste0("R²=", round(r2, 2),
                     "\nRMSE=", round(rmse, 3)),
      hjust = 1.1, vjust = -0.5
    ) +
    ggplot2::theme_classic()
}

save_plot <- function(plot_obj, filename, w = 7, h = 5) {
  ggplot2::ggsave(filename, plot_obj, width = w, height = h)
}

plot_variance <- function(r) {
  terra::plot(r, col = viridisLite::viridis(256))
}

plot_reference_vs_prediction <- function(ref, pred) {

  df <- as.data.frame(c(ref, pred), xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("x", "ref", "pred")

  ggplot2::ggplot(df, ggplot2::aes(ref, pred)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_abline(slope = 1, linetype = "dashed") +
    ggplot2::theme_classic()
}
