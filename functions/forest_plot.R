# forest plot created by ggplot2 ----

#' @param data `list` a list of data frames
#' @param lag `char` lag time window
#' @param title `char` title of the plot
#' @param subtitle `char` subtitle of the plot
geom_forest <- function(data, lag, title, subtitle) {
  data <- lapply(data, function(x) subset(x, lag_time == lag))

  df_os <- data[["OS"]] %>%
    transform(
      p_f = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)),
      index = seq(nrow(.), 1)
    ) %>%
    {
      inset(
        .,
        which(extract(., "lower.95") == 0),
        c("HR", "lower.95", "upper.95"),
        c(NA, NA, NA)
      )
    } %>%
    {
      inset(
        .,
        which(is.na(extract(., "lower.95"))),
        c("HR_f", "p_f"),
        c("-", "-")
      )
    }

  df_css <- data[["CSS"]] %>%
    transform(
      p_f = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)),
      index = seq(nrow(.), 1)
    ) %>%
    {
      inset(
        .,
        which(extract(., "lower.95") == 0),
        c("HR", "lower.95", "upper.95"),
        c(NA, NA, NA)
      )
    } %>%
    {
      inset(
        .,
        which(is.na(extract(., "lower.95"))),
        c("HR_f", "p_f"),
        c("-", "-")
      )
    }

  x_lmt <- c(
    min(df_os$lower.95, df_css$lower.95, na.rm = TRUE),
    max(df_os$HR, df_css$HR, na.rm = TRUE) + 0.1
  )

  p1 <- ggplot(data = df_os) +
    geom_text(
      aes(
        x = 0,
        y = index,
        label = unlist(cancer_names[cancer_site])
      ),
      hjust = 0,
      na.rm = TRUE
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = "Cancer type") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0, face = "bold"))

  p2 <- ggplot(data = df_os) +
    geom_text(
      aes(
        x = 0,
        y = index,
        label = HR_f
      ),
      hjust = 0.5,
      na.rm = TRUE
    ) +
    labs(title = "HR (95% CI)\nfor OS") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin()
    )

  p3 <- ggplot(data = df_os) +
    geom_text(
      aes(
        x = 0,
        y = index,
        label = p_f
      ),
      hjust = 0.5,
      na.rm = TRUE
    ) +
    labs(title = expression(bolditalic("p"))) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin()
    )

  p4 <- ggplot(data = df_os) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(
      aes(
        x = HR,
        y = index,
        xmin = lower.95,
        xmax = upper.95
      ),
      size = 0.1,
      na.rm = TRUE
    ) +
    geom_segment(
      aes(
        x = x_lmt[2],
        y = c(
          index[which(upper.95 > x_lmt[2])],
          rep(NA, times = length(index) - length(which(upper.95 > x_lmt[2])))
        ),
        xend = x_lmt[2] + 0.1,
        yend = c(
          index[which(upper.95 > x_lmt[2])],
          rep(NA, times = length(index) - length(which(upper.95 > x_lmt[2])))
        )
      ),
      arrow = arrow(length = unit(3, "points"), type = "closed"),
      na.rm = TRUE
    ) +
    xlab("Hazard ratio (95% CI)") +
    coord_cartesian(xlim = x_lmt) +
    scale_x_continuous(expand = c(0, 0.1)) +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(r = 10)
    )

  p5 <- ggplot(data = df_css) +
    geom_text(
      aes(
        x = 0,
        y = index,
        label = HR_f
      ),
      hjust = 0.5,
      na.rm = TRUE
    ) +
    labs(title = "HR (95% CI)\nfor CSS") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin()
    )

  p6 <- ggplot(data = df_css) +
    geom_text(
      aes(
        x = 0,
        y = index,
        label = p_f
      ),
      hjust = 0.5,
      na.rm = TRUE
    ) +
    labs(title = expression(bolditalic("p"))) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin()
    )

  p7 <- ggplot(data = df_css) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    geom_pointrange(
      aes(
        x = HR,
        y = index,
        xmin = lower.95,
        xmax = upper.95
      ),
      size = 0.1,
      na.rm = TRUE
    ) +
    geom_segment(
      aes(
        x = x_lmt[2],
        y = c(
          index[which(upper.95 > x_lmt[2])],
          rep(NA, times = length(index) - length(which(upper.95 > x_lmt[2])))
        ),
        xend = x_lmt[2] + 0.1,
        yend = c(
          index[which(upper.95 > x_lmt[2])],
          rep(NA, times = length(index) - length(which(upper.95 > x_lmt[2])))
        )
      ),
      arrow = arrow(length = unit(3, "points"), type = "closed"),
      na.rm = TRUE
    ) +
    xlab("Hazard ratio (95% CI)") +
    coord_cartesian(xlim = x_lmt) +
    scale_x_continuous(expand = c(0, 0.1)) +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(r = 10)
    )

  p_all <- p1 + p2 + p3 + p4 + p5 + p6 + p7 +
    plot_layout(nrow = 1, width = c(2, 1.2, 0.8, 2, 1.2, 0.8, 2)) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    )

  p_all
}
