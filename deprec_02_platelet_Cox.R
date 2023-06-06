
# env settings ------------------------------------------------------------

library(magrittr)
library(survival)
library(openxlsx)
dir.create("02", FALSE)
load("01/tidy_data_diagnosis_after_attending.RData")
load("00/ICD_of_cancers.RData")
load("00/ICD_of_cancers_excl_bood.RData")
load("00/fun_UKb_Cox_loop.RData")
load("00/fun_UKb_Cox_loop_excluded.RData")


# Cox regression for platelet ---------------------------------------------

tidy_data_dia_after_att <- tidy_data_dia_after_att %>% 
  transform(
    PlateletPer10 = Platelet / 10,
    PlateletPer100 = Platelet / 100
  )

## platelet (continuous per 1, per 10, per 100, 300 cutoff, and 400 cutoff)
## lag time (0 days, 6 months i.e. 365.25/2 days, and 1 year i.e. 365.25 days)
platelet_multi_Cox <- UKb_Cox_loop(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = c(
    "Age_at_recruitment",
    "Sex",
    "asprin",
    "SMOKING_STATUS",
    "ALCOHOL_STATUS",
    "BMI",
    "Townsend_deprivation_index.TDI._at_recruitment",
    "Race"
  ),
  ICD = "ICD",
  ICDrexp = ICD_rexp,
  attending = "Date_attending",
  lag = c(0, 365.25 / 2, 365.25),
  target = c("Platelet", "PlateletPer10", "PlateletPer100", "Platelet300", "Platelet400")
)

platelet_uni_Cox <- UKb_Cox_loop(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = NULL,
  ICD = "ICD",
  ICDrexp = ICD_rexp,
  attending = "Date_attending",
  lag = c(0, 365.25 / 2, 365.25),
  target = c("Platelet", "PlateletPer10", "PlateletPer100", "Platelet300", "Platelet400")
)

platelet_multi_Cox_excl_blood <- UKb_Cox_loop_excluded(
  data = tidy_data_dia_after_att,
  start = "Date_of_diagnosis",
  end = "Date_end",
  event = "OS",
  eventcode = 1,
  covars = c(
    "Age_at_recruitment", 
    "Sex", 
    "asprin", 
    "SMOKING_STATUS", 
    "ALCOHOL_STATUS", 
    "BMI", 
    "Townsend_deprivation_index.TDI._at_recruitment",
    "Race"
  ),
  ICD = "ICD",
  ICDrexp = ICD_rexp_excl_bood,
  attending = "Date_attending",
  lag = c(0, 365.25/2, 365.25),
  target = c("Platelet", "PlateletPer10", "PlateletPer100", "Platelet300", "Platelet400")
)

platelet_uni_Cox_excl_blood <- UKb_Cox_loop_excluded(
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
  target = c("Platelet", "PlateletPer10", "PlateletPer100", "Platelet300", "Platelet400")
)

# save data ---------------------------------------------------------------

write.xlsx(platelet_multi_Cox, "02/platelet_multi_Cox.xlsx", TRUE)
write.xlsx(platelet_uni_Cox, "02/platelet_uni_Cox.xlsx", TRUE)
write.xlsx(platelet_multi_Cox_excl_blood, "02/platelet_multi_Cox_excl_blood.xlsx", TRUE)
write.xlsx(platelet_uni_Cox_excl_blood, "02/platelet_uni_Cox_excl_blood.xlsx", TRUE)