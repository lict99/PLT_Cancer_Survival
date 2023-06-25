# env settings ----

library(survival)
library(openxlsx)
library(magrittr)

load("00/cancer_ICD_codes_with_attr.RData")
load("01/whole_cancer_data_for_Cox.RData")
source("functions/Cox_regression.R")

dir.create("06", FALSE)

# Cox regression for aspirin ----

multi_lag <- list(
  c(-Inf, 0),
  c(0, Inf),
  c(365.25 / 2, Inf),
  c(365.25, Inf)
)

# vars_per_100 <- c(
#   "platelet_per_100",
#   "age",
#   "sex",
#   "aspirin",
#   "smoking_status",
#   "alcohol_status",
#   "bmi",
#   "TDI",
#   "race",
#   "fu_time",
#   "fu_event"
# )

# vars_300 <- vars_per_100 %>% inset(1, "platelet_300")
# vars_400 <- vars_per_100 %>% inset(1, "platelet_400")

# ## multivariate
# platelet_per_100_multi <- lapply(
#   whole_cancer_data,
#   function(x) {
#     run_Cox_loop(
#       data_list = x,
#       mlagtime = multi_lag,
#       vars = vars_per_100,
#       target = "platelet"
#     )
#   }
# )

## univariate
aspirin_uni <- lapply(
  whole_cancer_data,
  function(x) {
    run_Cox_loop(
      data_list = x,
      mlagtime = multi_lag,
      vars = c("aspirin", "fu_time", "fu_event"),
      target = "aspirin"
    )
  }
)

# # data saving ----

# for (i in ls(pattern = "platelet")) {
#   write.xlsx(
#     get(i), paste0("02/", i, ".xlsx"), TRUE
#   )
# }
