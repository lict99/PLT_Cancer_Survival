# env settings ----

library(magrittr)
library(survival)
library(openxlsx)

load("00/cancer_ICD_codes_with_attr.RData")
load("01/whole_cancer_data_for_Cox.RData")
source("functions/Cox_regression.R")

dir.create("02", FALSE)

# Cox regression ----

multi_lag <- list(
  c(0, Inf),
  c(365.25 / 2, Inf),
  c(365.25, Inf)
)

vars_per_100 <- c(
  "platelet_per_100",
  "age",
  "sex",
  "asprin",
  "smoking_status",
  "alcohol_status",
  "bmi",
  "TDI",
  "race",
  "fu_time",
  "fu_event"
)

vars_300 <- vars_per_100 %>% inset(1, "platelet_300")
vars_400 <- vars_per_100 %>% inset(1, "platelet_400")

## multivariate
platelet_per_100_multi <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = vars_per_100,
      target = "platelet"
    )
  }
)

platelet300_multi <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = vars_300,
      target = "platelet"
    )
  }
)

platelet400_multi <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = vars_400,
      target = "platelet"
    )
  }
)

## univariate
platelet_per_100_uni <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("platelet_per_100", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

platelet300_uni <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("platelet_300", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

platelet400_uni <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("platelet_400", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

# data saving ----

for (i in ls(pattern = "platelet")) {
  write.xlsx(
    get(i), paste0("02/", i, ".xlsx"), TRUE
  )
}
