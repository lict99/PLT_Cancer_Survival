# env settings -----------------------------------------------------------------

library(magrittr)

load("01/UKB_all_info.RData")
load("00/cancer_ICD_codes_with_attr.RData")

source("functions/different_survival.R", local = TRUE)
source("functions/Cox_regression.R", local = TRUE)

# expiration of research -------------------------------------------------------

date_end <- as.Date("2021-07-01")

# extract different survival outcomes ------------------------------------------

UKb_diagnosis_cancer <- subset(
  UKb_diagnosis,
  grepl(cancer_ICD_codes[["All_sites"]], ICD_diagnosis)
)

whole_cancer_data_pre <- list(
  OS = extract_survival(
    death_codes = lapply(cancer_ICD_codes, function(x) "."),
    diag_df = UKb_diagnosis_cancer
  ),
  CSS = extract_survival(
    death_codes = cancer_ICD_codes,
    diag_df = UKb_diagnosis_cancer
  )
)

## include cancer types which have 1000 cases or more
## only consider eligible cases
## (cancer diagnosis after assessment canter visit)
n1k <- lapply(
  whole_cancer_data_pre,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = "eid",
      lagtime = c(0, Inf)
    ) %>%
      lapply(nrow) %>%
      unlist() %>%
      extract(is_weakly_greater_than(., 1000)) %>%
      sort(decreasing = TRUE) %>%
      names()
  }
)

whole_cancer_data <- mapply(
  function(x, y) {
    x[y]
  },
  whole_cancer_data_pre,
  n1k,
  SIMPLIFY = FALSE
)

# data saving ------------------------------------------------------------------

save(
  whole_cancer_data,
  file = "01/whole_cancer_data_for_Cox.RData",
  compress = FALSE
)
save(n1k, file = "01/cancers_with_more_than_1k_cases.RData", compress = FALSE)
