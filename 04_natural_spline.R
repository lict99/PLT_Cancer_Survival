# env settings ----

library(magrittr)
library(survival)
library(Hmisc)
library(smoothHR)
library(lmtest)
library(ggplot2)
library(scales)
library(patchwork)

load("01/whole_cancer_data_for_Cox.RData")
load("00/cancer_ICD_codes_with_attr.RData")
load("00/cancer_names.RData")
source("functions/Cox_regression.R")
source("functions/natural_spline.R")

dir.create("04", FALSE)

# calculation and plotting ----

lag <- c(0, Inf)

vars <- c(
  "fu_time", "fu_event",
  "platelet", "age",
  "sex",
  "aspirin", "smoking_status", "alcohol_status", "bmi", "TDI", "race"
)

plot_list <- list()
for (i in names(whole_cancer_data)) {
  plot_list[[i]] <- lapply(
    extract_Cox_data(
      data_list = whole_cancer_data[[i]],
      lagtime = lag,
      vars = vars
    ),
    function(x) try(cal_hr(x), silent = TRUE)
  ) %>%
    lapply(function(x) try(geom_hr(x), silent = TRUE)) %>%
    extract(., which(sapply(., is.list))) %>%
    {
      list <- list()
      for (j in names(.)) {
        list[[j]] <- extract2(., j) +
          labs(
            title = cancer_names[[j]],
            caption = paste(
              "The reference value is the median of platelet counts",
              "\n",
              paste("Cancer survival:", i),
              sep = ""
            )
          )
      }
      list
    }
}

p_patch <- plot_list[["OS"]][["All_sites"]] +
  plot_list[["CSS"]][["All_sites"]] +
  guide_area() +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    design = c(
      area(1, 1, 9, 5),
      area(1, 6, 9, 10),
      area(10, 1, 10, 10)
    ),
    guides = "collect"
  ) &
  theme(legend.position = "bottom")

# plots saving ----

for (i in names(plot_list)) {
  for (j in names(plot_list[[i]])) {
    ggsave(
      paste0("04/", i, "_", cancer_names[[j]], ".pdf"),
      plot_list[[i]][[j]],
      height = 9,
      width = 12,
      device = grDevices::cairo_pdf
    )
  }
}

ggsave(
  "04/Pan-cancer.pdf",
  p_patch,
  height = 9,
  width = 18,
  device = grDevices::cairo_pdf
)
