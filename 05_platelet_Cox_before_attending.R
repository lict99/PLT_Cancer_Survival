# env settings ----

library(magrittr)
library(survival)
library(openxlsx)

load("00/cancer_ICD_codes_with_attr.RData")
load("01/whole_cancer_data_for_Cox.RData")
source("functions/Cox_regression.R")

dir.create("05", FALSE)

# Cox regression for diagnosis before attending ----

lag <- list(c(-Inf, 0))

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

## model 2
platelet_per_100_bam2 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = lag,
      vars = vars_per_100,
      target = "platelet"
    )
  }
)

platelet300_bam2 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = lag,
      vars = vars_300,
      target = "platelet"
    )
  }
)

platelet400_bam2 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = lag,
      vars = vars_400,
      target = "platelet"
    )
  }
)

## model 1
platelet_per_100_bam1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = lag,
      vars = c("platelet_per_100", "age", "sex", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

platelet300_bam1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = lag,
      vars = c("platelet_300", "age", "sex", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

platelet400_bam1 <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = lag,
      vars = c("platelet_400", "age", "sex", "fu_time", "fu_event"),
      target = "platelet"
    )
  }
)

# data saving ----

for (i in ls(pattern = "^platelet.+_bam")) {
  write.xlsx(
    get(i), paste0("05/", i, ".xlsx"), TRUE
  )
}
