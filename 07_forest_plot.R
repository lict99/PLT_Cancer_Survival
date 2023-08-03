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
    "10\u00b9\u00b9/L"
  } else if (grepl("300_", i)) {
    "\u2265300\u00d710\u2079/L"
  } else if (grepl("400_", i)) {
    "\u2265400\u00d710\u2079/L"
  } else {
    stop("Invalid objects! Check ls() call.")
  }
  model <- switch(strsplit(i, "_")[[1]][length(strsplit(i, "_")[[1]])],
    "m1" = "model 1",
    "m2" = "model 2",
    stop("Invalid model name!")
  )
  for (j in lag_time) {
    fp <- geom_forest(
      get(i),
      j,
      paste0(
        "Cox proportional hazards model",
        " (", model, ") ",
        "\nfor platelet counts ",
        " (", pc_type, ") ",
        "on cancer survival"
      )
    )
    ggsave(
      paste0(
        "07/",
        i,
        "_",
        strsplit(j, " ")[[1]][1],
        ".pdf"
      ),
      fp,
      width = 13,
      height = 9,
      units = "in",
      device = grDevices::cairo_pdf
    )
  }
}
