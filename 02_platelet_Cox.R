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

multi_var1 <- c(
  "platelet_per_100",
  "age",
  "sex",
  "asprin",
  "smoking_status",
  "alcohol_status",
  "bmi",
  "TDI",
  "race"
)

multi_var2 <- multi_var1 %>% inset(1, "platelet_300")
multi_var3 <- multi_var1 %>% inset(1, "platelet_400")

## multivariate
platelet_per_100_multi_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = multi_var1
)

platelet300_multi_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = multi_var2
)

platelet400_multi_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = multi_var3
)

## univariate
platelet_per_100_uni_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = "platelet_per_100"
)

platelet300_uni_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = "platelet_300"
)

platelet400_uni_Cox <- run_Cox_loop(
  mlagtime = multi_lag,
  vars = "platelet_400"
)

# data arrangement --------------------------------------------------------


