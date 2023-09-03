# env settings ----

library(magrittr)
library(survival)
library(openxlsx)

load("00/cancer_ICD_codes_with_attr.RData")
load("01/whole_cancer_data_for_Cox.RData")
source("functions/Cox_regression.R")

dir.create("03", FALSE)

# Cox regression ----

## multiple lag time (0 day, 6 months, and 1 year)
multi_lag <- list(
  c(0, Inf),
  c(365.25 / 2, Inf),
  c(365.25, Inf)
)

## variables to be analyzed in complex Cox regression
vars_per_100 <- c(
  "platelet_per_100",
  "age",
  "sex",
  "aspirin",
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

## complex model
platelet_per_100_m2 <- lapply(
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

platelet300_m2 <- lapply(
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

platelet400_m2 <- lapply(
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

## basic model
platelet_per_100_m1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("platelet_per_100", "age", "sex", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

platelet300_m1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("platelet_300", "age", "sex", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

platelet400_m1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("platelet_400", "age", "sex", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

# data saving ----

save(list = ls(pattern = "^platelet.+_m"), file = "03/platelet_Cox.RData")

for (i in ls(pattern = "^platelet.+_m")) {
  write.xlsx(
    get(i), paste0("03/", i, ".xlsx"), TRUE
  )
}
