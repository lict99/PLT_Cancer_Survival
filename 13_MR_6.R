# env settings ----

library(magrittr)
library(TwoSampleMR)
library(openxlsx)

load("09/MR_harmonised_data.RData")
load("09/snp_proxies.RData")
snp6 <- read.xlsx("src/SNP/PC_SNPs_6.xlsx", 1) %>%
  use_series(snp) %>%
  as.character()

dir.create("13", FALSE)

# MR ----

snp_need_prx <- which(!is.element(snp6, mr_data[[1]][[1]]$SNP)) %>%
  extract(snp6, .)

is.element(c("affy20084012", "affy28520978"), c(names(ukb_snp_ind), names(ukb_snp_ind_supp)))
