# env settings ------------------------------------------------------------

library(magrittr)

load("01/UKB_all_info.RData")
load("00/cancer_ICD_codes_with_attr.RData")

# end of follow-up time ---------------------------------------------------
date_end <- as.Date("2021-07-01")

# cancer diagnosis --------------------------------------------------------

UKb_diagnosis_cancer <- subset(
  UKb_diagnosis,
  grepl(cancer_ICD_codes[["All_sites"]], ICD_diagnosis)
)

diag_per_cancer <- lapply(
  cancer_ICD_codes,
  function(x) {
    if (attr(x, "incl") == TRUE) {
      subset(UKb_diagnosis_cancer, grepl(x, ICD_diagnosis))
    } else if (attr(x, "incl") == FALSE) {
      subset(UKb_diagnosis_cancer, !grepl(x, ICD_diagnosis))
    } else {
      NA
    }
  }
)

# cancer death ------------------------------------------------------------

UKb_death_cancer <- subset(
  UKb_death,
  grepl(cancer_ICD_codes[["All_sites"]], ICD10_death)
)

death_per_cancer <- lapply(
  cancer_ICD_codes,
  function(x) {
    if (attr(x, "incl") == TRUE) {
      subset(UKb_death_cancer, grepl(x, ICD10_death))
    } else if (attr(x, "incl") == FALSE) {
      subset(UKb_death_cancer, !grepl(x, ICD10_death))
    } else {
      NA
    }
  }
)

# index of diagnosis to death ---------------------------------------------

death_fu_end <- UKb_death[, -3] %>%
  extract(order(use_series(., eid), use_series(., date_death)), ) %>%
  extract(!duplicated(use_series(., eid)), )

cancer_diag_to_death <- list()
for (i in names(cancer_ICD_codes)) {
  diag_i <- diag_per_cancer[[i]] %>%
    extract(order(use_series(., eid), use_series(., date_diagnosis)), ) %>%
    extract(!duplicated(use_series(., eid)), -2) %>%
    set_colnames(c("eid", "date_cancer_diagnosis"))

  death_i <- death_per_cancer[[i]] %>%
    extract(order(use_series(., eid), use_series(., date_death)), ) %>%
    extract(!duplicated(use_series(., eid)), -3) %>%
    set_colnames(c("eid", "date_cancer_death"))

  cancer_diag_to_death[[i]] <- merge(
    diag_i, death_i,
    by = "eid", all.x = TRUE
  ) %>%
    merge(death_fu_end, by = "eid", all.x = TRUE) %>%
    transform(
      cancer_death = ifelse(is.na(date_cancer_death), 0, 1),
      fu_end = ifelse(
        is.na(date_death),
        as.character(date_end),
        as.character(date_death)
      ) %>%
        as.Date()
    )
}

cancer_data_Cox <- lapply(
  cancer_diag_to_death,
  function(x) {
    merge(x, UKb_baseline, by = "eid", all.x = TRUE)
  }
)

whole_cancer_data <- cancer_data_Cox %>%
  lapply(
    function(x) {
      transform(
        x,
        lag_time = difftime(
          date_cancer_diagnosis,
          date_attending,
          units = "days"
        ),
        fu_time = difftime(
          fu_end,
          date_cancer_diagnosis,
          units = "days"
        ),
        platelet_per_100 = platelet / 100,
        platelet_300 = ifelse(platelet >= 300, "YES", "NO") %>% factor(),
        platelet_400 = ifelse(platelet >= 400, "YES", "NO") %>% factor()
      )
    }
  )

# data saving -------------------------------------------------------------

save(whole_cancer_data, file = "01/whole_cancer_data_for_Cox.RData")
