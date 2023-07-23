# PRS plot by ggplot2 ----

gg_prs <- function(
    data,
    annotate = NULL,
    x_col,
    y_col = "platelet",
    covars = c("age", "sex")) {
  ## calculation
  cor <- cor(
    data[, x_col],
    data[, y_col],
    method = "pearson"
  ) %>%
    sprintf("%.3f", .)

  cor_p <- cor.test(
    data[, x_col],
    data[, y_col],
    method = "pearson",
    alternative = "two.sided"
  ) %>%
    extract2("p.value") %>%
    {
      ifelse(is_less_than(., 0.001), sprintf("%.3e", .), sprintf("%.3f", .))
    }


  fit_smr <- lm(
    as.formula(
      paste(
        paste(y_col, "~", x_col),
        paste(covars, collapse = "+"),
        sep = "+"
      )
    ),
    data
  ) %>%
    summary()

  f <- fit_smr$fstatistic[1] %>% sprintf("%.3f", .)
  fv1 <- fit_smr$fstatistic[2]
  fv2 <- fit_smr$fstatistic[3]
  fit_p <- pf(
    fit_smr$fstatistic[1],
    fit_smr$fstatistic[2],
    fit_smr$fstatistic[3],
    lower.tail = FALSE
  ) %>%
    {
      ifelse(is_less_than(., 0.001), sprintf("%.3e", .), sprintf("%.3f", .))
    }

  ## plotting
  col_d <- "#7876B1FF"
  col_l <- "#BC3C29FF"

  p_x <- ggplot() +
    geom_density(
      aes(x = data[, x_col]),
      fill = alpha(col_d, 0.8),
      color = col_d
    ) +
    labs(
      x = NULL,
      y = if (x_col == "prs_w") {
        "Density of Weighted PRS"
      } else if (x_col == "prs_u") {
        "Density of Unweighted PRS"
      } else {
        stop("x_col should be either \'prs_w\' or \'prs_u\'!")
      }
    ) +
    theme_classic() +
    theme(
      plot.margin = margin(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.length.x = unit(0, "pt")
    )

  p_y <- ggplot() +
    geom_density(
      aes(x = data[, y_col]),
      fill = alpha(col_d, 0.8),
      color = col_d
    ) +
    labs(x = NULL, y = "Density of Platelet Count") +
    coord_flip() +
    theme_classic() +
    theme(
      plot.margin = margin(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.length.y = unit(0, "pt")
    )

  p_xy <- ggplot() +
    geom_point(
      aes(x = data[, x_col], y = data[, y_col]),
      alpha = 0.15,
      size = 0.1
    ) +
    geom_smooth(
      aes(x = data[, x_col], y = data[, y_col]),
      method = "lm",
      formula = "y ~ x",
      color = col_l
    ) +
    xlab(
      if (x_col == "prs_w") {
        "Weighted PRS"
      } else if (x_col == "prs_u") {
        "Unweighted PRS"
      } else {
        stop("x_col should be either \'prs_w\' or \'prs_u\'!")
      }
    ) +
    ylab(expression(paste("Platelet Count (", 10^9, "/L)", sep = ""))) +
    annotate(
      "text",
      x = max(data[, x_col]),
      y = max(data$platelet),
      label =
        paste(
          annotate,
          "\n",
          "Sample size = ", nrow(data),
          "\n",
          "Pearson r = ", cor, "; P-value = ", cor_p,
          "\n",
          "F(", fv1, ", ", fv2, ") = ", f, "; P-value = ", fit_p,
          sep = ""
        ),
      hjust = 1,
      vjust = 1
    ) +
    theme_classic()

  p_all <- p_x + plot_spacer() + p_xy + p_y +
    plot_layout(
      ncol = 2,
      byrow = TRUE,
      widths = c(8, 2),
      heights = c(2, 8)
    )
  p_all
}
