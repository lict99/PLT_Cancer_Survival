# env settings ----

rm(list = ls())
gc()

library(magrittr)
library(survival)
library(openxlsx)
library(data.table)

load("01/whole_cancer_data_for_Cox.RData")
load("00/cancer_ICD_codes_with_attr.RData")
source("functions/Cox_regression.R")

ukb_snp_ind <- fread("src/SNP/UKb_snp_individual.csv", data.table = FALSE)
ukb_snp_ind_supp <- fread(
  "src/SNP/UKb_snp_individual_supp.csv",
  data.table = FALSE
)
ukb_snp_sum <- fread("src/SNP/UKb_snp_summary.txt", data.table = FALSE)

pc_snp <- read.xlsx("src/SNP/PC_SNPs.xlsx", 1)
pc_snp_prx <- read.xlsx("src/SNP/query_snp.xlsx", 1)

dir.create("08", FALSE)

# data preprocessing ----

## transform SNP summary information
snp_smr <- ukb_snp_sum %>%
  transform(
    affy = paste0("affy", affy_id),
    snp = paste0("rs", rs_id)
  ) %>%
  subset(subset = is.element(snp, union(pc_snp[, 1], pc_snp_prx[, 1]))) %>%
  subset(
    subset = is.element(
      affy,
      union(colnames(ukb_snp_ind), colnames(ukb_snp_ind_supp))
    )
  ) %>%
  transform(
    ref_allele = apply(
      ., 1,
      function(x) {
        if (x["ref"] == "A") {
          x["allele_a"]
        } else if (x["ref"] == "B") {
          x["allele_b"]
        } else {
          NA
        }
      }
    ),
    eff_allele = apply(
      ., 1,
      function(x) {
        if (x["ref"] == "A") {
          x["allele_b"]
        } else if (x["ref"] == "B") {
          x["allele_a"]
        } else {
          NA
        }
      }
    )
  )

## extract baseline data for Cox regression
Cox_data <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = c("eid", "fu_time", "fu_event", "age", "sex"),
      lagtime = c(0, Inf)
    )
  }
)

## extract all patients with cancer
## and corresponding SNP information
cancer_eid <- Cox_data[["OS"]][["All_sites"]][, "eid"]

snp_ind_ca <- subset(
  merge(ukb_snp_ind, ukb_snp_ind_supp, by = "eid", all = TRUE),
  is.element(eid, cancer_eid)
)

## substitute genotype with 0, 1, 2
## with a safe process of substitution
snp_ind <- data.frame(eid = snp_ind_ca$eid)
for (i in snp_smr$affy) {
  ref_i <- snp_smr[snp_smr$affy == i, "ref_allele"]
  eff_i <- snp_smr[snp_smr$affy == i, "eff_allele"]
  homo_ref_i <- paste(ref_i, ref_i)
  hete_i <- c(paste(ref_i, eff_i), paste(eff_i, ref_i))
  homo_eff_i <- paste(eff_i, eff_i)
  col_snp <- snp_ind_ca[, i]
  for (j in seq_along(col_snp)) {
    if (!is.na(col_snp[j])) {
      if (col_snp[j] == homo_ref_i) {
        col_snp[j] <- 0
      } else if (col_snp[j] == hete_i[1]) {
        col_snp[j] <- 1
      } else if (col_snp[j] == hete_i[2]) {
        col_snp[j] <- 1
      } else if (col_snp[j] == homo_eff_i) {
        col_snp[j] <- 2
      } else {
        col_snp[j] <- NA
      }
    }
  }
  snp_ind[, i] <- as.numeric(col_snp)
}

## merge baseline data and SNP data
data_Cox_snp <- list()
for (i in names(Cox_data)) {
  data_Cox_snp[[i]] <- lapply(
    Cox_data[[i]],
    function(x) {
      merge(x, snp_ind, by = "eid", all.x = TRUE)
    }
  )
}

# Cox regression of SNP ----

## compute the association between SNPs and cancer survival
## using UKB individual data
Cox_snp <- list()
for (j in names(data_Cox_snp)) {
  Cox_snp[[j]] <- lapply(
    data_Cox_snp[[j]],
    function(x) {
      snp_df <- data.frame()
      for (i in snp_smr$affy) {
        if (is.element("sex", colnames(x))) {
          fml <- as.formula(
            paste("Surv(fu_time, fu_event == 1) ~ age + sex +", i)
          )
          fit <- coxph(
            fml,
            data = x,
            singular.ok = TRUE,
            x = TRUE
          )
        } else {
          fml <- as.formula(
            paste("Surv(fu_time, fu_event == 1) ~ age +", i)
          )
          fit <- coxph(
            fml,
            data = x,
            singular.ok = TRUE,
            x = TRUE
          )
        }
        smr <- summary(fit)
        df <- data.frame(
          affy = i,
          snp = snp_smr[snp_smr$affy == i, "snp"],
          ref = snp_smr[snp_smr$affy == i, "ref_allele"],
          eff = snp_smr[snp_smr$affy == i, "eff_allele"],
          coef = smr[["coefficients"]][i, "coef"],
          se = smr[["coefficients"]][i, "se(coef)"],
          p = smr[["coefficients"]][i, "Pr(>|z|)"],
          n = fit[["n"]],
          eaf = (sum(fit[["x"]][, i])) / (fit[["n"]] * 2)
        )
        snp_df <- rbind(snp_df, df)
      }
      snp_df
    }
  )
}

# data saving ----

save(snp_smr, file = "08/pc_snp_summary.RData")
save(snp_ind_ca, file = "08/genotype_with_cancer.RData")
save(Cox_snp, file = "08/Cox_regression_of_SNP.RData")
