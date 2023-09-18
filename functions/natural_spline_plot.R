# HR estimates ----

#' @param data `data.frame` data to run Cox proportional hazards model
cal_hr <- function(data) {
  if (any(grepl("sex", colnames(data)))) {
    fit <- coxph(
      formula = Surv(fu_time, fu_event == 1) ~
        rcspline.eval(platelet, nk = 3, inclx = TRUE) + sex +
        age + aspirin + smoking_status + alcohol_status + bmi + TDI + race,
      data = data,
      x = TRUE,
      singular.ok = TRUE
    )
  } else {
    fit <- coxph(
      formula = Surv(fu_time, fu_event == 1) ~
        rcspline.eval(platelet, nk = 3, inclx = TRUE) +
        age + aspirin + smoking_status + alcohol_status + bmi + TDI + race,
      data = data,
      x = TRUE,
      singular.ok = TRUE
    )
  }

  fit2 <- coxph(
    Surv(fu_time, fu_event == 1) ~ .,
    data = data,
    x = TRUE,
    singular.ok = TRUE
  )

  hr <- smoothHR(data = data, coxfit = fit)
  hr_point <- predict.HR(
    hr,
    predictor = "platelet",
    prob = 0.5,
    prediction.values = seq(
      min(hr$dataset$platelet),
      max(hr$dataset$platelet),
      length.out = 100
    ),
    conf.level = 0.95
  ) %>%
    as.data.frame() %>%
    within({
      HR <- exp(LnHR)
      lower_95 <- exp(`lower .95`)
      upper_95 <- exp(`upper .95`)
    })

  list(
    dataset = hr$dataset,
    sumr = summary(fit),
    hr_point = hr_point,
    lrtest = lrtest(fit2, fit)[2, "Pr(>Chisq)"]
  )
}

# natural spline plot created by ggplot2 ----

#' @param data `list` data from cal_hr()
geom_hr <- function(data) {
  dataset <- data$dataset
  p <- data$lrtest
  point <- data$hr_point
  dst <- density(dataset$platelet)
  dst_lmt <- c(0, max(dst$y))
  hr_lmt <- c(0, max(point$HR))
  dst_y_rsc <- rescale(dst$y, hr_lmt, dst_lmt)
  col_dst <- "#7876B1FF"
  col_hr <- "#BC3C29FF"
  col_ci <- "black"

  p1 <- ggplot() +
    geom_ribbon(
      aes(
        x = dst$x,
        ymin = min(dst_y_rsc),
        ymax = dst_y_rsc,
        fill = "Density"
      ),
      color = col_dst
    ) +
    geom_hline(yintercept = 1, color = "gray") +
    geom_line(
      aes(x = platelet, y = lower_95, color = "95% confidence interval"),
      point,
      linetype = 2
    ) +
    geom_line(
      aes(x = platelet, y = upper_95, color = "95% confidence interval"),
      point,
      linetype = 2
    ) +
    geom_line(
      aes(x = platelet, y = HR, color = "Hazard ratio"),
      point,
      linewidth = 0.9,
      linetype = 1
    ) +
    annotate(
      "label",
      x = min(dst$x),
      y = ceiling(hr_lmt[2]),
      label = ifelse(
        p < 0.001,
        paste(
          "paste(italic(P), '-', 'value')==",
          paste0(
            "'", unlist(strsplit(sprintf("%.3e", p), "e"))[1], "'",
            "%*%10^", unlist(strsplit(sprintf("%.3e", p), "e"))[2]
          ),
          sep = ""
        ),
        paste(
          "paste(italic(P), '-', 'value')==",
          "'", sprintf("%.3f", p), "'",
          sep = ""
        )
      ),
      label.size = 0,
      parse = TRUE,
      hjust = 0,
      vjust = 0.1
    ) +
    scale_x_continuous(
      name = expression(paste("Platelet counts (", 10^9, "/L)"))
    ) +
    scale_y_continuous(
      breaks = seq(0, ceiling(hr_lmt[2]), by = 1),
      name = "Hazard ratio",
      sec.axis = sec_axis(
        trans = ~ rescale(., dst_lmt, c(0, max(.))),
        name = "Density of platelet counts"
      )
    ) +
    scale_color_manual(
      breaks = c("Hazard ratio", "95% confidence interval"),
      values = c(col_hr, col_ci),
      guide = guide_legend(
        title = NULL,
        override.aes = list(linetype = c(1, 2), linewidth = 0.6)
      )
    ) +
    scale_fill_manual(
      name = NULL,
      breaks = c("Density"),
      values = c(alpha(col_dst, 0.8))
    ) +
    coord_cartesian(
      ylim = c(0, ceiling(hr_lmt[2]))
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black")
    )

  p1
}

# multiple natural spline plots into one plot ----

#' @param plots `list` a list of plots
#' @param caption `char` caption of the plot
geom_multi_hr <- function(
    plots,
    caption) {
  p <- wrap_plots(
    lapply(
      plots,
      function(x) x + theme(plot.caption = element_blank())
    )
  ) +
    guide_area() +
    plot_annotation(
      tag_levels = NULL,
      caption = caption,
      theme = theme(plot.caption = element_text(size = 15))
    ) +
    plot_layout(
      ncol = 4,
      guides = "collect"
    ) &
    theme(
      legend.position = "right",
      legend.text = element_text(size = 15)
    )

  p
}
