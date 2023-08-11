# env settings ----

library(magrittr)
library(ggplot2)
library(patchwork)

load("00/cancer_names.RData")
load("03/platelet_Cox.RData")
source("functions/forest_plot.R")

dir.create("07", FALSE)

# forest plot ----

lag_time <- c(
  "0 to Inf day(s)", "182.625 to Inf day(s)", "365.25 to Inf day(s)"
)

for (i in ls(pattern = "platelet.+m[12]")) {
  pc_type <- if (grepl("_100_", i)) {
    "100\u00d710\u2079/L"
  } else if (grepl("300_", i)) {
    "\u2265300\u00d710\u2079/L"
  } else if (grepl("400_", i)) {
    "\u2265400\u00d710\u2079/L"
  } else {
    stop("Invalid objects! Check ls() call.")
  }
  model <- switch(strsplit(i, "_")[[1]][length(strsplit(i, "_")[[1]])],
    "m1" = "Adjusted for age and sex",
    "m2" = paste(
      "Adjusted for",
      "age, sex,",
      "race, BMI, TDI, aspirin use, smoking status, and alcohol status"
    ),
    stop("Invalid model name!")
  )
  for (j in lag_time) {
    lagtime <- switch(strsplit(j, " ")[[1]][1],
      `365.25` = " with a lag time of 1 year",
      `182.625` = " with a lag time of 6 months",
      `0` = "",
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
      paste0("07/", i, "_", strsplit(j, " ")[[1]][1], ".pdf"),
      fp,
      width = 13,
      height = 9,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }
}
