# ----

arrange_lag <- function(data) {
  vars <- lapply(data, names) %>% unlist()
  if (length(unique(vars)) != length(vars) / length(data)) {
    stop("Column names are not identical!")
  }
  plot_data <- mapply(
    function(x, nm) {
      transform(x, survival = nm)
    },
    data,
    names(data),
    SIMPLIFY = FALSE
  )
  Reduce(rbind, plot_data)
}

tmp <- arrange_lag(platelet_per_100_m2) %>%
  transform(
    cancer_site = sapply(cancer_site, function(x) cancer_names[[x]]) %>%
      factor(levels = unique(.)),
    lag_time = sapply(
      lag_time,
      function(x) as.numeric(strsplit(x, "-")[[1]][[1]]) / 365.25
    ),
    survival = factor(
      survival,
      levels = c("OS", "CSS"),
      labels = c("Overall survival", "Cancer-specific survival")
    )
  )
# ----
geom_lag_trend <- function(data) {
  p1 <- ggplot(data) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray") +
    geom_linerange(
      aes(x = lag_time, ymin = lower.95, ymax = upper.95, color = survival),
      alpha = 0.5
    ) +
    geom_line(aes(x = lag_time, y = HR, color = survival), alpha = 0.5) +
    scale_x_continuous(
      name = "Lag time (year)",
      breaks = c(0, 0.5, 1, 3, 5),
      labels = c("≥0", "≥0.5", "≥1", "≥3", "≥5")
    ) +
    facet_wrap(vars(cancer_site), scales = "free_y") +
    theme(
      panel.grid = element_blank(),
      legend.position = "top"
    )
  p1
}

geom_lag_trend(tmp)
