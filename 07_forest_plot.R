# env settings ----

library(magrittr)
library(ggplot2)
library(patchwork)

load("00/cancer_names.RData")
load("03/platelet_Cox.RData")
source("functions/forest_plot.R", local = TRUE)

dir.create("07", FALSE)

# forest plot ----

## multiple lag time
lag_time <- c(
  "0-Inf", "182.625-Inf", "365.25-Inf", "1095.75-Inf"
)

## visualize the results of Cox regression by forest plot
## there are some display problems using Cairo graphics device
## but it is not a big deal
for (i in ls(pattern = "platelet.+m[12]")) {
  pc_type <- if (grepl("_100_", i)) {
    "100\u00d710\u2079/L increase"
  } else if (grepl("300_", i)) {
    "\u2265300\u00d710\u2079/L vs. \u003c300\u00d710\u2079/L"
  } else if (grepl("400_", i)) {
    "\u2265400\u00d710\u2079/L vs. \u003c400\u00d710\u2079/L"
  } else {
    stop("Invalid objects! Check ls() call.")
  }
  model <- switch(strsplit(i, "_")[[1]][length(strsplit(i, "_")[[1]])],
    m1 = "Adjusted for age and sex",
    m2 = paste(
      "Adjusted for",
      "age, sex,",
      "race, BMI, TDI, aspirin use, smoking status, and alcohol status"
    ),
    stop("Invalid model name!")
  )
  for (j in lag_time) {
    lagtime <- switch(j,
      `0-Inf` = "",
      `182.625-Inf` = " with a lag time of ≥6 months",
      `365.25-Inf` = " with a lag time of ≥1 year",
      `1095.75-Inf` = " with a lag time of ≥3 years",
      stop("Invalid lag time!")
    )
    fp <- geom_forest(
      data = get(i),
      lag = j,
      title = paste0(
        "Cox proportional hazards model",
        "\nfor platelet counts ",
        " (", pc_type, ") ",
        "on cancer survival"
      ),
      subtitle = paste0(model, lagtime)
    )
    ggsave(
      paste0("07/", i, "_", j, ".pdf"),
      fp,
      width = 11.5,
      height = 7.5,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }
}
