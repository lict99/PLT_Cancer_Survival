# env settings ------------------------------------------------------------

library(survival)
library(magrittr)

load("01/whole_cancer_data_for_Cox.RData")
load("00/cancer_ICD_codes_with_attr.RData")
load("00/functions.RData")

UKb_snp_ind <- read.csv("src/SNP/UKb_snp_individual.csv")
UKb_snp_sum <- read.delim("src/SNP/UKb_snp_summary.txt")
snp_platelet <- read.csv("src/SNP/snp_platelet.csv")

dir.create("09", FALSE)

# data preprocessing ------------------------------------------------------

snp_sum <- UKb_snp_sum %>%
  transform(
    affy = paste0("affy", affy_id),
    snp = paste0("rs", rs_id)
) %>%
  subset(subset = is.element(affy, colnames(UKb_snp_ind))) %>%
  transform(
    allele_ref = apply(
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
    allele_eff = apply(
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

cancer_data_extr <- extract_Cox_data(
  data_list = whole_cancer_data,
  vars = c("eid", "fu_time", "cancer_death", "age", "sex"),
  lagtime = 0
)

## a safe process of substitution
UKb_snp_ind_ca <- subset(
  UKb_snp_ind, 
  is.element(eid, cancer_data_extr$All_sites$eid)
)
snp_ind <- data.frame(eid = UKb_snp_ind_ca$eid)
for (i in snp_sum$affy) {
  ref_i <- snp_sum[snp_sum$affy == i, "allele_ref"]
  eff_i <- snp_sum[snp_sum$affy == i, "allele_eff"]
  homo_ref_i <- paste(ref_i, ref_i)
  hete_i <- c(paste(ref_i, eff_i), paste(eff_i, ref_i))
  homo_eff_i <- paste(eff_i, eff_i)
  col_snp <- UKb_snp_ind_ca[, i]
  col_val <- numeric()
  for (j in seq_len(length(col_snp))) {
    if (col_snp[j] == homo_ref_i) {
      col_val[j] <- 0
    } else if (col_snp[j] == hete_i[1]) {
      col_val[j] <- 1
    } else if (col_snp[j] == hete_i[2]) {
      col_val[j] <- 1
    } else if (col_snp[j] == homo_eff_i) {
      col_val[j] <- 2
    } else {
      col_val[j] <- NA
    }
  }
  snp_ind[, i] <- col_val
}

data_Cox_snp <- lapply(
  cancer_data_extr,
  function(x) {
    merge(x, snp_ind, by = "eid", all.x = TRUE)
  }
)

# Cox regression of SNP ---------------------------------------------------

Cox_snp <- list()
for (i in names(cancer_ICD_codes)) {
  data <- data_Cox_snp[[i]]
  
  y <- Surv(
    time = data[, "fu_time"], 
    event = data[, "cancer_death"]
  )
  
  snp_df <- data.frame()
  for (j in snp_sum$affy) {
    if (all(is.na(data[, j]))) {
      next 
    }
    
    if (length(attr(cancer_ICD_codes[[i]], "sex")) == 0) {
      fit <- coxph(
        as.formula(paste("y ~ age + sex", j, sep = " + ")),
        data = data,
        singular.ok = TRUE,
        model = TRUE,
        x = TRUE
      )
    } else if (length(attr(cancer_ICD_codes[[i]], "sex")) == 1) {
      fit <- coxph(
        as.formula(paste("y ~ age", j, sep = " + ")),
        data = data,
        singular.ok = TRUE,
        model = TRUE,
        x = TRUE
      )
    }
    
    summ <- summary(fit)
    
    df <- data.frame(
      affy = j,
      snp = snp_sum[snp_sum$affy == j, "snp"],
      ref = snp_sum[snp_sum$affy == j, "allele_ref"],
      eff = snp_sum[snp_sum$affy == j, "allele_eff"],
      coef = summ[["coefficients"]][j, "coef"],
      se = summ[["coefficients"]][j, "se(coef)"],
      p = summ[["coefficients"]][j, "Pr(>|z|)"],
      maf = (sum(fit[["x"]][, j])) / (fit[["n"]] * 2),
      ph_pass = tryCatch(
        ifelse(all(cox.zph(fit)[["table"]][, "p"] > 0.05), "YES", "NO"), 
        error = function(e) {NA_character_}
      )
    )
    
    snp_df <- rbind(snp_df, df)
  }
  
  Cox_snp[[i]] <- snp_df
}

# data saving -------------------------------------------------------------

save(snp_ind, file = "09/SNP_individual_with_numeric_values.RData")
save(Cox_snp, file = "09/Cox_regression_of_SNP.RData")
