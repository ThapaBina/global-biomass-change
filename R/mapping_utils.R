# =========================================================
# MAPPING UTILITIES
# =========================================================

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

