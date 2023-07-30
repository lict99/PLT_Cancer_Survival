# env settings ----

library(magrittr)
library(forestplot)

load("00/cancer_names.RData")
load("03/platelet_Cox.RData")

dir.create("07", FALSE)

# data preparation ----

## extract and format data
prep_fp_data <- function(data_df,
                         colnames = c(
                           "cancer_site", "HR_f", "p_f",
                           "HR", "lower.95", "upper.95"
                         ),
                         lag,
                         full_names = cancer_names) {
  df <- subset(data_df, lag_time == lag) %>%
    extract(, colnames, drop = FALSE)
  if (identical(df$cancer_site, names(full_names))) {
    df$cancer_site <- unlist(full_names)
  }
  df
}
pdf("07/forest_plot.pdf", width = 13, height = 13 * 0.618)
(function() {
  tmp <- prep_fp_data(platelet300_m2[[1]], lag = "0 to Inf day(s)")
  tmp2 <- prep_fp_data(platelet300_m2[[2]], lag = "0 to Inf day(s)")

  p1 <- forestplot(
    x = tmp,
    mean = HR,
    lower = lower.95,
    upper = upper.95,
    labeltext = c(cancer_site, HR_f, p_f),
    zero = 1,
    boxsize = 0.1,
    clip = c(floor(min(tmp$HR)), ceiling(max(tmp$HR)))
  ) %>%
    fp_add_header(
      cancer_site = "Cancer type",
      HR_f = "HR (95% CI) for OS",
      p_f = "P-value"
    )

  p2 <- forestplot(
    x = tmp2,
    mean = HR,
    lower = lower.95,
    upper = upper.95,
    labeltext = c(HR_f, p_f),
    zero = 1,
    boxsize = 0.1,
    align = c("r", "r"),
    clip = c(floor(min(tmp2$HR)), ceiling(max(tmp2$HR)))
  ) %>%
    fp_add_header(
      HR_f = "HR (95% CI) for CSS",
      p_f = "P-value"
    )
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = c(4.3, 3))))
  pushViewport(viewport(layout.pos.col = 1))
  plot(p1)
  popViewport(1)
  pushViewport(viewport(layout.pos.col = 2))
  plot(p2)
  popViewport(2)
})()
dev.off()
