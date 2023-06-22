# env settings ----
library(magrittr)

load("01/UKB_all_info.RData")
load("00/cancer_ICD_codes_with_attr.RData")
source("functions/different_survival.R")

# end of follow-up time ----
date_end <- as.Date("2021-07-01")

# extract different survival outcomes ----

whole_cancer_data <- list(
  OS = extract_survival(
    death_codes = lapply(cancer_ICD_codes, function(x) ".")
  ),
  CSS = extract_survival(
    death_codes = cancer_ICD_codes
  )
  # death_by_blood_and_circulation = extract_survival(
  #   death_codes = lapply(
  #     cancer_ICD_codes,
  #     function(x) {
  #       "^D5[0-35-9]|^D[68][0-9]|^D7[0-7]|^I|^28[0-9]|^39[0-9]|^4[0-5][0-9]"
  #     }
  #   )
  # )
)

# data saving ----

save(whole_cancer_data, file = "01/whole_cancer_data_for_Cox.RData")
