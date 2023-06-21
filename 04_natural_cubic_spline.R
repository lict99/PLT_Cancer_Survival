# env settings ----

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

# calculation and plotting ----
lag <- c(0, Inf)

vars <- c(
  "fu_time", "fu_event",
  "platelet", "age",
  "sex",
  "asprin", "smoking_status", "alcohol_status", "bmi", "TDI", "race"
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
              "The reference value is the median of platelet counts.",
              "\n",
              paste("survival:", i),
              "; ",
              paste("lag time:", paste(lag, collapse = " to "), "day(s)"),
              sep = ""
            )
          )
      }
      list
    }
}
