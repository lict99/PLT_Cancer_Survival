gg_forest <- function(data, lag, title) {
  data <- lapply(data, function(x) subset(x, lag_time == lag))
  p1 <- ggplot(data = data[["OS"]]) +
    geom_text(
      aes(
        x = 0,
        y = seq(nrow(data[["OS"]]), 1),
        label = unlist(cancer_names[cancer_site])
      ),
      hjust = 0
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = "Cancer type") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.1, face = "bold"))

  p2 <- ggplot(data = data[["OS"]]) +
    geom_text(
      aes(
        x = 0,
        y = seq(nrow(data[["OS"]]), 1),
        label = HR_f
      ),
      hjust = 0.5
    ) +
    labs(title = "HR (95% CI)\nfor OS") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p3 <- ggplot(data = data[["OS"]]) +
    geom_text(
      aes(
        x = 0,
        y = seq(nrow(data[["OS"]]), 1),
        label = p_f
      ),
      hjust = 0.5
    ) +
    labs(title = "P-value") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p4 <- ggplot(data = data[["OS"]]) +
    geom_pointrange(
      aes(
        x = HR,
        y = seq(nrow(data[["OS"]]), 1),
        xmin = lower.95,
        xmax = upper.95
      ),
      size = 0.1
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    xlab("Hazard ratio (95% CI)") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin()
    )

  p5 <- ggplot(data = data[["CSS"]]) +
    geom_text(
      aes(
        x = 0,
        y = seq(nrow(data[["CSS"]]), 1),
        label = HR_f
      ),
      hjust = 0.5
    ) +
    labs(title = "HR (95% CI)\nfor CSS") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p6 <- ggplot(data = data[["CSS"]]) +
    geom_text(
      aes(
        x = 0,
        y = seq(nrow(data[["CSS"]]), 1),
        label = p_f
      ),
      hjust = 0.5
    ) +
    labs(title = "P-value") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p7 <- ggplot(data = data[["CSS"]]) +
    geom_pointrange(
      aes(
        x = HR,
        y = seq(nrow(data[["CSS"]]), 1),
        xmin = lower.95,
        xmax = upper.95
      ),
      size = 0.1
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    xlab("Hazard ratio (95% CI)") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin()
    )

  p_all <- p1 + p2 + p3 + p4 + p5 + p6 + p7 +
    plot_layout(nrow = 1, width = c(1.8, 1, 1, 2, 1, 1, 2)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    )
  p_all
}
