# env settings ----

library(survival)
library(openxlsx)
library(magrittr)

load("00/cancer_ICD_codes_with_attr.RData")
load("01/whole_cancer_data_for_Cox.RData")
source("functions/Cox_regression.R")

dir.create("06", FALSE)

# Cox regression for aspirin ----

## multiple lag time
multi_lag <- list(
  c(-Inf, 0),
  c(0, Inf),
  c(365.25 / 2, Inf),
  c(365.25, Inf)
)

## variables to be analyzed in complex Cox regression
vars_aspirin <- c(
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

## basic model
## extract the effect of aspirin on cancer survival
aspirin_m1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("aspirin", "age", "sex", "fu_time", "fu_event"),
      target = "aspirin"
    )
  }
)

## complex model
## extract the effect of aspirin on cancer survival
aspirin_m2 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = vars_aspirin,
      target = "aspirin"
    )
  }
)

# data saving ----

for (i in ls(pattern = "^aspirin_m")) {
  write.xlsx(
    get(i), paste0("06/", i, ".xlsx"), TRUE
  )
}
