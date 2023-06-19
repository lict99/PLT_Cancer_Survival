# env settings ------------------------------------------------------------

library(cutoff)
library(survival)
library(survminer)
load("01/tidy_data_diagnosis_after_attending.RData")


# cutoff calculation ------------------------------------------------------

tidy_data_dia_after_att <- transform(
  tidy_data_dia_after_att,
  futime = as.numeric(Date_end - Date_of_diagnosis)
)

data1 <- tidy_data_dia_after_att[, c("futime", "OS", "Platelet")] %>%
  na.omit()

cutoff <- cox(
  data = data1,
  time = "futime",
  y = "OS",
  x = "Platelet",
  cut.numb = 1,
  n.per = 0.25,
  y.per = 0.10
)

cutoff2 <- surv_cutpoint(
  data1,
  "futime",
  "OS",
  "Platelet"
)
