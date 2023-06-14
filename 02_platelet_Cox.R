# env settings ------------------------------------------------------------

library(magrittr)
library(survival)
library(openxlsx)

load("00/cancer_ICD_codes_with_attr.RData")
load("01/whole_cancer_data_for_Cox.RData")
# load("00/functions.RData")

dir.create("02", FALSE)

# Cox regression ----------------------------------------------------------

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
  "cancer_death"
)

vars_300 <- vars_per_100 %>% inset(1, "platelet_300")
vars_400 <- vars_per_100 %>% inset(1, "platelet_400")

## multivariate
platelet_per_100_multi_Cox <- extract_Cox_data(
  vars = vars_per_100,
  lagtime = multi_lag[[1]]
) %>%
  run_Cox_per_lag() %>%
  extract_smr_data()

platelet300_multi_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = vars_300,
  target = "platelet"
)

platelet400_multi_Cox

## univariate
platelet_per_100_uni_Cox

platelet300_uni_Cox

platelet400_uni_Cox

# data arrangement --------------------------------------------------------
