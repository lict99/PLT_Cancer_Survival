# env settings ------------------------------------------------------------

library(magrittr)

dir.create("01", FALSE)

# data input --------------------------------------------------------------

csv_files <- list.files(path = "src", pattern = "\\.csv")
names_csv <- gsub("\\.csv", "", csv_files)
for (i in seq_along(csv_files)) {
  assign(
    names_csv[i],
    read.csv(file.path("src", csv_files[i]))
  )
}

# data arrangement ---------------------------------------------------------

## basic info
Basic2 <- subset(
  Basic_information,
  subset = !duplicated(eid),
  select = -Smoking_status
)
rm(Basic_information)
gc()

## cancer diagnosis
Date_cancer2 <- subset(
  Date_cancer_registry,
  !is.na(Date_of_cancer_diagnosis)
) %>%
  transform(
    ICD_diagnosis = paste0(
      use_series(., Type_of_cancer_ICD10),
      use_series(., Type_of_cancer_ICD9)
    ) %>%
      gsub("NA", "", .)
  )
rm(Date_cancer_registry)
gc()

diagnosis_ICD10_2 <- subset(
  diagnosis_ICD10,
  (DIAGNOSES_ICD10 != "") & (DATE_ICD10 != "")
)
rm(diagnosis_ICD10)
gc()

diagnosis_ICD9_2 <- subset(
  diagnosis_ICD9,
  (DIAGNOSES_ICD9 != "") & (DATE_ICD9 != "")
) %>%
  transform(
    DIAGNOSES_ICD9 = ifelse(
      grepl("^[A-Z]", DIAGNOSES_ICD9, ignore.case = TRUE),
      paste(DIAGNOSES_ICD9, "icd9", sep = "-"),
      DIAGNOSES_ICD9
    )
  )
rm(diagnosis_ICD9)
gc()

## combination of data regarding diagnosis
ICD_diagnoses <- rbind(
  Date_cancer2[, c(1, 5, 2)] %>% set_colnames(c("eid", "diagnosis", "date_of_diagnosis")),
  diagnosis_ICD10_2 %>% set_colnames(c("eid", "diagnosis", "date_of_diagnosis")),
  diagnosis_ICD9_2 %>% set_colnames(c("eid", "diagnosis", "date_of_diagnosis"))
) %>%
  extract(
    order(
      use_series(., eid),
      use_series(., date_of_diagnosis)
    ),
  ) %>%
  subset(!duplicated(.)) %>%
  set_rownames(NULL)
rm(Date_cancer2, diagnosis_ICD10_2, diagnosis_ICD9_2)
gc()

## weight
weight2 <- subset(weight, !duplicated(eid))
rm(weight)
gc()

# mergence of all data ----------------------------------------------------

data_merged <- merge(
  Basic2,
  asprin,
  by = "eid",
  all = TRUE
) %>%
  merge(
    Date_attending,
    by = "eid",
    all = TRUE
  ) %>%
  merge(
    Platelet,
    by = "eid",
    all = TRUE
  ) %>%
  merge(
    smoking_alcohol_day,
    by = "eid",
    all = TRUE
  ) %>%
  merge(
    weight2[, -2],
    by = "eid",
    all = TRUE
  ) %>%
  merge(
    basic_disease,
    by = "eid",
    all = TRUE
  ) %>%
  set_colnames(
    tolower(colnames(.))
  )

# transformation of data --------------------------------------------------

tidy_colnames <- colnames(data_merged) %>%
  inset(c(2, 4), c("age", "TDI"))

UKb_baseline <- data_merged %>%
  set_colnames(tidy_colnames) %>%
  transform(
    sex = factor(sex, levels = c("Female", "Male")),
    race = race %>%
      ifelse(equals(., "Do not know"), NA, .) %>%
      ifelse(equals(., "Prefer not to answer"), NA, .),
    asprin = factor(asprin, levels = c(0, 1), labels = c("NO", "YES")),
    date_attending = as.Date(date_attending),
    smoking_status = factor(
      smoking_status,
      levels = c(0, 1, 2),
      labels = c("Never", "Previous", "Current")
    ),
    alcohol_status = factor(
      alcohol_status,
      levels = c(0, 1, 2),
      labels = c("Never", "Previous", "Current")
    ),
    diabetes = factor(diabetes, c(0, 1), c("NO", "YES")),
    dyslipidemia = factor(dyslipidemia, c(0, 1), c("NO", "YES")),
    hypertension = factor(hypertension, c(0, 1), c("NO", "YES"))
  ) %>%
  transform(
    race = ifelse(race == "British", "British", "Others") %>%
      factor(levels = c("Others", "British"))
  )

UKb_diagnosis <- ICD_diagnoses %>%
  set_colnames(c("eid", "ICD_diagnosis", "date_diagnosis")) %>%
  transform(date_diagnosis = as.Date(date_diagnosis))

UKb_death <- subset(
  Date_death,
  (Date_death != "") & (ICD10 != "")
) %>%
  set_colnames(c("eid", "date_death", "ICD10_death")) %>%
  transform(date_death = as.Date(date_death))

# data saving -------------------------------------------------------------

save(UKb_baseline, UKb_diagnosis, UKb_death, file = "01/UKB_all_info.RData")
