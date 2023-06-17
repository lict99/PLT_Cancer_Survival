# env settings ------------------------------------------------------------

library(magrittr)
library(survival)
library(Hmisc)
library(smoothHR)
library(lmtest)
library(ggplot2)
library(ggsci)
library(scales)

load("01/whole_cancer_data_for_Cox.RData")
load("00/cancer_ICD_codes_with_attr.RData")
source("functions/Cox_regression.R")
source("functions/natural_cubic_spline.R")

dir.create("04", FALSE)

# calculation -------------------------------------------------------------

vars <- c(
  "fu_time", "cancer_death",
  "platelet", "age",
  "sex",
  "asprin", "smoking_status", "alcohol_status", "bmi", "TDI", "race"
)

data <- extract_Cox_data(lagtime = c(0, Inf), vars = vars)

hr_smth <- lapply(data, function(x) try(cal_hr(x), silent = TRUE))

# plotting ----------------------------------------------------------------

plot_list <- lapply(hr_smth, function(x) try(geom_hr(x), silent = TRUE)) %>%
  extract(., which(sapply(., is.list)))

# mplot <- plot_grid(
#   plotlist = plot_list,
#   align = "none",
#   ncol = sqrt(length(plot_list)) %>% ceiling(),
#   labels = names(plot_list) %>%
#     gsub("_", " ", ., fixed = TRUE) %>%
#     gsub("E ", "Excluding ", ., fixed = TRUE),
#   label_fontface = "plain",
#   hjust = 0,
#   vjust = 1.1,
#   label_x = 0.1
# )

# save_plot("04/rcs_plot_lag0.pdf", mplot, base_height = 30)
