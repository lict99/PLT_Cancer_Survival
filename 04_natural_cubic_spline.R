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
load("00/cancer_names.RData")
source("functions/Cox_regression.R")
source("functions/natural_cubic_spline.R")

dir.create("04", FALSE)

# calculation -------------------------------------------------------------

vars <- c(
  "fu_time", "fu_event",
  "platelet", "age",
  "sex",
  "asprin", "smoking_status", "alcohol_status", "bmi", "TDI", "race"
)

hr_smth <- lapply(
  extract_Cox_data(
    data_list = whole_cancer_data[["OS"]],
    lagtime = c(0, Inf),
    vars = vars
  ),
  function(x) try(cal_hr(x), silent = TRUE)
)

# plotting ----------------------------------------------------------------

plot_list <- lapply(hr_smth, function(x) try(geom_hr(x), silent = TRUE)) %>%
  extract(., which(sapply(., is.list))) %>%
  {
    list <- list()
    for (i in names(.)) {
      list[[i]] <- extract2(., i) +
        labs(
          title = cancer_names[[i]],
          caption = "The reference value is the median of platelet counts"
        )
    }
    list
  }
