
# env settings ------------------------------------------------------------

library(survival)
library(openxlsx)
library(magrittr)
load("01/tidy_data_diagnosis_after_attending.RData")
load("00/fun_UKb_Cox_loop_excluded.RData")
load("00/fun_UKb_Cox_loop.RData")
load("00/ICD_of_cancers_excl_bood.RData")
load("00/ICD_of_cancers.RData")
dir.create("06", FALSE)

# Cox regression for asprin -----------------------------------------------

bsc_dses <- read.csv("src/basic_disease.csv")
tidy_data_dia_after_att <- merge(
  tidy_data_dia_after_att,
  bsc_dses,
  by = "eid",
  all.x = TRUE,
  sort = FALSE
) %>% 
  transform(
    Diabetes = factor(Diabetes, levels = c(0, 1), labels = c("NO", "YES")),
    Dyslipidemia = factor(Dyslipidemia, levels = c(0, 1), labels = c("NO", "YES")),
    Hypertension = factor(Hypertension, levels = c(0, 1), labels = c("NO", "YES"))
  )

## asprin
## lag time (0 days, 6 months i.e. 365.25/2 days, and 1 year i.e. 365.25 days)
asprin_multi_Cox <- UKb_Cox_loop(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = c(
    "Age_at_recruitment", 
    "Sex", 
    "Platelet", 
    "SMOKING_STATUS", 
    "ALCOHOL_STATUS", 
    "BMI", 
    "Townsend_deprivation_index.TDI._at_recruitment",
    "Race",
    "Diabetes",
    "Dyslipidemia",
    "Hypertension"
  ),
  ICD = "ICD",
  ICDrexp = ICD_rexp,
  attending = "Date_attending",
  lag = c(0, 365.25/2, 365.25),
  target = c("asprin")
)

asprin_uni_Cox <- UKb_Cox_loop(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = NULL,
  ICD = "ICD",
  ICDrexp = ICD_rexp,
  attending = "Date_attending",
  lag = c(0, 365.25/2, 365.25),
  target = c("asprin")
)

asprin_multi_Cox_excl_blood <- UKb_Cox_loop_excluded(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = c(
    "Age_at_recruitment", 
    "Sex", 
    "Platelet", 
    "SMOKING_STATUS", 
    "ALCOHOL_STATUS", 
    "BMI", 
    "Townsend_deprivation_index.TDI._at_recruitment",
    "Race",
    "Diabetes",
    "Dyslipidemia",
    "Hypertension"
  ),
  ICD = "ICD",
  ICDrexp = ICD_rexp_excl_bood,
  attending = "Date_attending",
  lag = c(0, 365.25/2, 365.25),
  target = c("asprin")
)

asprin_uni_Cox_excl_blood <- UKb_Cox_loop_excluded(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = NULL,
  ICD = "ICD",
  ICDrexp = ICD_rexp_excl_bood,
  attending = "Date_attending",
  lag = c(0, 365.25/2, 365.25),
  target = c("asprin")
)

# save data ---------------------------------------------------------------

write.xlsx(asprin_multi_Cox, "06/asprin_multi_Cox.xlsx", TRUE)
write.xlsx(asprin_multi_Cox_excl_blood, "06/asprin_multi_Cox_excl_blood.xlsx", TRUE)
write.xlsx(asprin_uni_Cox, "06/asprin_uni_Cox.xlsx", TRUE)
write.xlsx(asprin_uni_Cox_excl_blood, "06/asprin_uni_Cox_excl_blood.xlsx", TRUE)
