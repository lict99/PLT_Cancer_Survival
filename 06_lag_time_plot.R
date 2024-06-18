# env settings -----------------------------------------------------------------

library(ggplot2)

load("00/cancer_names.RData")
load("01/cancers_with_more_than_1k_cases.RData")
load("03/platelet_Cox.RData")

dir.create("06", FALSE)

for (data_nm in ls(pattern = "platelet.+m[12]")) {
  pc_type <- if (grepl("_100_", data_nm)) {
    "100\u00d710\u2079/L increase"
  } else if (grepl("300_", data_nm)) {
    "\u2265300\u00d710\u2079/L vs. \u003c300\u00d710\u2079/L"
  } else if (grepl("400_", data_nm)) {
    "\u2265400\u00d710\u2079/L vs. \u003c400\u00d710\u2079/L"
  } else {
    stop("Invalid objects! Check ls() call.")
  }
  model <- strsplit(data_nm, "_")[[1]][length(strsplit(data_nm, "_")[[1]])]
  model_lab <- switch(model,
    m1 = "Adjusted for age and sex",
    m2 = paste(
      "Adjusted for",
      "age, sex,",
      "ethnic background, BMI, TDI, smoking status, and alcohol status"
    ),
    stop("Invalid model name!")
  )
  data <- get(data_nm)
  for (surv in names(data)) {
    surv_lab <- switch(surv,
      OS = "overall survival",
      CSS = "cancer-specific survival"
    )
    d <- data[[surv]]
    d <- d[d$lag_time != "0-Inf", , drop = FALSE]
    d$lag_time <- factor(
      d$lag_time,
      levels = c(
        "0-182.625", "182.625-365.25", "365.25-1095.75", "1095.75-Inf"
      ),
      labels = c("0-6m", "6m-1y", "1-3y", "3+y")
    )
    d$cancer_site <- factor(
      d$cancer_site,
      levels = n1k$OS,
      labels = vapply(n1k$OS, function(x) cancer_names[[x]], character(1L))
    )
    p <- ggplot(d) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
      geom_errorbar(
        aes(x = lag_time, ymin = lower.95, ymax = upper.95),
        width = 0.2
      ) +
      geom_point(aes(x = lag_time, y = HR)) +
      facet_wrap(vars(cancer_site), scales = "free_y") +
      labs(
        title = paste0(
          "Cox proportional hazards model",
          "\nfor platelet counts ",
          " (", pc_type, ") ",
          "on ", surv_lab
        ),
        subtitle = paste0(model_lab),
        x = "Lag time interval",
        y = "Hazard ratio (95% CI)"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    ggsave(
      paste0("06/", data_nm, "_", surv, ".pdf"),
      p,
      device = grDevices::cairo_pdf,
      width = 11,
      height = 8
    )
  }
}
