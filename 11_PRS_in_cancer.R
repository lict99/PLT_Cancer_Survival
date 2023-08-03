# env settings ----

library(magrittr)
library(survival)
library(ggplot2)
library(patchwork)

load("01/whole_cancer_data_for_Cox.RData")
load("08/genotype_with_cancer.RData")
load("08/pc_snp_summary.RData")
load("10/MR_harmonised_data.RData")
load("00/cancer_ICD_codes_with_attr.RData")
load("00/cancer_names.RData")

source("functions/Cox_regression.R")
source("functions/PRS_plot.R")

dir.create("11", FALSE)

# PRS calculation ----

## extract IV information

iv_info <- mr_data[["OS"]][["All_sites"]] %>%
  subset(mr_keep == TRUE) %>%
  merge(
    snp_smr[, c("affy", "snp")],
    by.x = "SNP",
    by.y = "snp",
    all.x = TRUE
  )

## calculate SNP effect allele score
snp_score <- data.frame(eid = snp_ind_ca$eid) %>%
  (function(x) {
    for (i in iv_info$SNP) {
      affy <- snp_smr[snp_smr$snp == i, "affy"]
      col_geno <- snp_ind_ca[, affy]
      snp_info <- subset(iv_info, SNP == i)
      eff <- snp_info[, "effect_allele.exposure"]
      oth <- snp_info[, "other_allele.exposure"]
      homo_oth <- paste(oth, oth)
      homo_eff <- paste(eff, eff)
      heto <- c(paste(oth, eff), paste(eff, oth))
      for (j in seq_along(col_geno)) {
        if (!is.na(col_geno[j])) {
          if (col_geno[j] == homo_oth) {
            col_geno[j] <- 0
          } else if (col_geno[j] == heto[1]) {
            col_geno[j] <- 1
          } else if (col_geno[j] == heto[2]) {
            col_geno[j] <- 1
          } else if (col_geno[j] == homo_eff) {
            col_geno[j] <- 2
          } else {
            col_geno[j] <- NA
          }
        }
      }
      x[, i] <- as.numeric(col_geno)
    }
    return(x)
  })()

## un-weighted and weighted PRS for various cancers
prs <- data.frame(
  eid = snp_score$eid,
  prs_u = rowSums(snp_score[, -1]),
  prs_w = (function() {
    if (identical(names(snp_score[, -1]), iv_info$SNP)) {
      a <- as.matrix(snp_score[, -1])
      b <- matrix(iv_info$beta.exposure, ncol = 1)
    } else {
      stop("Not indentical in names!")
    }
    return(a %*% b)
  })()
)

prs_cancer <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = c("eid", "fu_time", "fu_event", "age", "sex", "platelet"),
      lagtime = c(0, Inf)
    ) %>%
      lapply(
        function(x) {
          merge(x, prs, by = "eid", all.x = TRUE) %>%
            na.omit()
        }
      )
  }
)

# plotting ----

prs_w_plot <- lapply(
  prs_cancer,
  function(x) {
    list <- list()
    for (i in names(x)) {
      list[[i]] <- geom_prs(
        data = x[[i]],
        annotate = cancer_names[[i]],
        x_col = "prs_w"
      )
    }
    list
  }
)

prs_u_plot <- lapply(
  prs_cancer,
  function(x) {
    list <- list()
    for (i in names(x)) {
      list[[i]] <- geom_prs(
        data = x[[i]],
        annotate = cancer_names[[i]],
        x_col = "prs_u"
      )
    }
    list
  }
)

# data saving ----

save(iv_info, file = "11/IV_info.RData")
