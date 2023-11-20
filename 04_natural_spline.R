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
source("functions/Cox_regression.R", local = TRUE)
source("functions/natural_spline_plot.R", local = TRUE)

dir.create("04", FALSE)

# calculation and plotting ----

## lag time of 0 day
## namely no lag time
lag <- c(0, Inf)

## variables to be analyzed in Cox regression with natural splines
vars <- c(
  "fu_time", "fu_event",
  "platelet", "age",
  "sex",
  "aspirin", "smoking_status", "alcohol_status", "bmi", "TDI", "race"
)

## all graphs in one list
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
              paste(
                "Cancer survival:",
                switch(i,
                  OS = "overall survival",
                  CSS = "cancer-specific survival"
                )
              ),
              sep = ""
            )
          )
      }
      list
    }
}

## combine pan-cancer plots
p_pan_cancer <- plot_list[["OS"]][["All_sites"]] +
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

## combine plots of overall survival
p_os <- geom_multi_hr(
  plots = plot_list[["OS"]][
    !is.element(
      names(plot_list[["OS"]]),
      c("All_sites")
    )
  ],
  caption = NULL
)

## combine plots of cancer-specific survival
p_css <- geom_multi_hr(
  plots = plot_list[["CSS"]][
    !is.element(
      names(plot_list[["CSS"]]),
      c("All_sites")
    )
  ],
  caption = NULL
)

# plots saving ----

ggsave(
  "04/Pan-cancer.pdf",
  p_pan_cancer,
  height = 6.5,
  width = 11,
  device = "pdf"
)

ggsave(
  "04/OS.pdf",
  p_os,
  height = 27,
  width = 27,
  device = "pdf"
)

ggsave(
  "04/CSS.pdf",
  p_css,
  height = 27,
  width = 27,
  device = "pdf"
)
