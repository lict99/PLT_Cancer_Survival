
# env settings ------------------------------------------------------------

library(survival)
library(survminer)
library(muhaz)
library(ggplot2)
library(cowplot)

load("01/whole_cancer_data_for_Cox.RData")
load("00/functions.RData")
load("00/cancer_ICD_codes_with_attr.RData")

dir.create("08", FALSE)

# calculation -------------------------------------------------------------


