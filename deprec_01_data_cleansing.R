
# env settings ------------------------------------------------------------

library(magrittr)
library(tidyverse)
dir.create("01", FALSE)
load("00/ICD_of_cancers.RData")
ref_dir <- "src/"

# read files --------------------------------------------------------------

csv_files <- list.files(path = ref_dir, pattern = "\\.csv")

names_csv <- gsub("\\.csv", "", csv_files)

for (i in seq_along(csv_files)) {
  assign(
    names_csv[i], 
    read.csv(paste0(ref_dir, csv_files[i]), check.names = T)
  )
}

# arrange data ------------------------------------------------------------

## basic info
BasicV2 <- subset(Basic_information, !duplicated(eid))

## cancer diagnosis
Date_cancerV2 <- subset(
  Date_cancer_registry, 
  !is.na(Date_of_cancer_diagnosis)
) %>% 
  transform(
    ICD = paste0(
      use_series(., Type_of_cancer_ICD10), 
      use_series(., Type_of_cancer_ICD9)
    ) %>% 
      gsub("NA", "", .)
  ) %>% 
  subset(
    grepl(ICD_rexp[[1]], ICD)
  )

diagnosis_ICD10V2 <- subset(diagnosis_ICD10, grepl(ICD_rexp[[1]], DIAGNOSES_ICD10))

diagnosis_ICD9V2 <- subset(diagnosis_ICD9, grepl(ICD_rexp[[1]], DIAGNOSES_ICD9))

## combine data of cancer diagnosis
## **仅保留样本第一次被诊断肿瘤的类型和日期**
cancer_info <- data.frame(
  eid = c(
    Date_cancerV2$eid,
    diagnosis_ICD10V2$eid,
    diagnosis_ICD9V2$eid
  ),
  ICD = c(
    Date_cancerV2$ICD,
    diagnosis_ICD10V2$DIAGNOSES_ICD10,
    diagnosis_ICD9V2$DIAGNOSES_ICD9
  ),
  Date_of_diagnosis = c(
    Date_cancerV2$Date_of_cancer_diagnosis,
    diagnosis_ICD10V2$DATE_ICD10,
    diagnosis_ICD9V2$DATE_ICD9
  )
) %>%
  magrittr::extract(
    order(
      use_series(., eid),
      use_series(., Date_of_diagnosis)
    ),
  ) %>% 
  subset(
    !duplicated(eid)
  )

## date of death
Date_deathV2 <- subset(Date_death, !duplicated(eid))

## weight
weightV2 <- subset(weight, !duplicated(eid))

# combine all data --------------------------------------------------------

tidy_data <- merge(
  cancer_info,
  BasicV2,
  by = "eid",
  all.x = T
) %>% 
  merge(
    asprin,
    by = "eid",
    all.x = T
  ) %>% 
  merge(
    Date_attending,
    by = "eid",
    all.x = T
  ) %>% 
  merge(
    Date_deathV2[,-3],
    by = "eid",
    all.x = T
  ) %>% 
  merge(
    Platelet,
    by = "eid",
    all.x = T
  ) %>% 
  merge(
    smoking_alcohol_day,
    by = "eid",
    all.x = T
  ) %>% 
  merge(
    weightV2[,-2],
    by = "eid",
    all.x = T
  ) %>% 
  merge(
    cancer_screening,
    by = "eid",
    all.x = T
  ) %>% 
  transform(
    Date_of_diagnosis = as.Date(Date_of_diagnosis),
    Age_at_recruitment65 = ifelse(
      Age_at_recruitment >= 65, "YES", "NO"
    ) %>% 
      as.factor(),
    Sex = factor(
      Sex,
      levels = c("Female", "Male")
    ),
    SMOKING_STATUS = factor(
      SMOKING_STATUS,
      levels = c(0, 1, 2),
      labels = c("Never", "Previous", "Current")
    ),
    ALCOHOL_STATUS = factor(
      ALCOHOL_STATUS, 
      levels = c(0, 1, 2), 
      labels = c("Never", "Previous", "Current")
    ),
    Race = na_if(Race, "Do not know") %>% 
      na_if("Prefer not to answer"),
    asprin = factor(
      asprin, 
      levels = c(0, 1), 
      labels = c("NO", "YES")
    ),
    Date_attending = as.Date(Date_attending),
    Date_end = ifelse(
      Date_death == "",
      "2021-07-01",
      Date_death
    ) %>% 
      as.Date(),
    OS = ifelse(
      Date_death == "", 0, 1
    ),
    Platelet300 = ifelse(
      Platelet >= 300, "YES" , "NO"
    ) %>% 
      as.factor(),
    Platelet400 = ifelse(
      Platelet >= 400, "YES", "NO"
    ) %>% 
      as.factor()
  ) %>% 
  transform(
    Race = ifelse(
      Race == "British", "British", "Others"
    ) %>% 
      as.factor()
  ) 

tidy_data_dia_after_att <- subset(
  tidy_data,
  Date_of_diagnosis >= Date_attending
)

# save data ---------------------------------------------------------------

save(tidy_data, file = "01/tidy_data.RData")
save(tidy_data_dia_after_att, file = "01/tidy_data_diagnosis_after_attending.RData")