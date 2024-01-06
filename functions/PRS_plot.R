# PRS plot created by ggplot2 ----

#' @param data `data.frame` data to visualize
#' @param annotate `char` annotation of the plot
#' @param x_col `char` column name of x-axis
#' @param y_col `char` column name of y-axis
geom_prs <- function(
    data,
    annotate = NULL,
    x_col,
    y_col = "platelet") {
  ## calculation
  cor <- cor(
    data[[x_col]],
    data[[y_col]],
    method = "pearson"
  ) %>%
    sprintf("%.3f", .)

  cor_p <- cor.test(
    data[[x_col]],
    data[[y_col]],
    method = "pearson",
    alternative = "two.sided"
  ) %>%
    extract2("p.value") %>%
    {
      ifelse(is_less_than(., 0.001), "< 0.001", paste("=", sprintf("%.3f", .)))
    }

  fit_smr <- lm(
    as.formula(paste(y_col, "~", x_col)),
    data
  ) %>%
    summary()

  r_squared <- fit_smr$r.squared %>% sprintf("%.3f", .)
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
      ifelse(is_less_than(., 0.001), "< 0.001", paste("=", sprintf("%.3f", .)))
    }

  ## plotting
  col_d <- "#7876B1FF"
  col_l <- "#BC3C29FF"

  plot_df <- data.frame(
    x = data[[x_col]],
    y = data[[y_col]]
  )

  p_x <- ggplot(data = plot_df) +
    geom_density(
      aes(x = x),
      fill = alpha(col_d, 0.8),
      color = col_d
    ) +
    labs(
      x = NULL,
      y = if (x_col == "prs_w") {
        "Density of weighted PRS"
      } else if (x_col == "prs_u") {
        "Density of unweighted PRS"
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

  p_y <- ggplot(data = plot_df) +
    geom_density(
      aes(x = y),
      fill = alpha(col_d, 0.8),
      color = col_d
    ) +
    labs(x = NULL, y = "Density of platelet counts") +
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

  p_xy <- ggplot(data = plot_df, aes(x = x, y = y)) +
    geom_hex(bins = 100) +
    scale_fill_gradient(low = "white", high = "black") +
    stat_ellipse(linetype = 2) +
    geom_smooth(
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
    ylab(expression(paste("Platelet counts (", 10^9, "/L)"))) +
    theme_classic() +
    theme(
      legend.position = c(0.9, 0.8),
      legend.justification = c(1, 0.5)
    )

  ann_df <- data.frame(
    x = rep(1, times = 5),
    y = 5:1,
    label = c(
      paste0("'", annotate, "'"),
      paste0("'Sample size'==", nrow(data)),
      paste0(
        "paste('Pearson r'=='", cor, "'",
        ",';'~~italic('p')*' ", cor_p, "')"
      ),
      paste0("'R'^2=='", r_squared, "'"),
      paste0(
        "paste(F[list(", fv1, ",", fv2, ")]=='", f, "'",
        ",';'~~italic('p')*' ", fit_p, "')"
      )
    )
  )

  p_ann <- ggplot(data = ann_df) +
    geom_text(
      aes(x = x, y = y, label = label),
      parse = TRUE,
      hjust = 0
    ) +
    coord_cartesian(xlim = c(1, 2)) +
    theme_void()

  p_xy_ann <- p_xy + inset_element(p_ann, 0.1, 0.65, 0.8, 0.9)

  p_all <- p_x + plot_spacer() + p_xy_ann + p_y +
    plot_layout(
      ncol = 2,
      byrow = TRUE,
      widths = c(8, 2),
      heights = c(2, 8),
      guides = "keep"
    )

  p_all
}
