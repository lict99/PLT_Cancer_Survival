
# env settings ------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(ggsci)
library(scales)
library(cowplot)

load("11/IV_info.RData")
load("11/functions_gg_prs.RData")
load("01/UKB_all_info.RData")
load("01/whole_cancer_data_for_Cox.RData")

ukb_snp_ind <- read.csv("src/SNP/UKb_snp_individual.csv")
ukb_snp_ind_supp <- read.csv("src/SNP/UKb_snp_individual_supp.csv")

dir.create("12", FALSE)

# calculation -------------------------------------------------------------

snp_ind <- merge(
  ukb_snp_ind, 
  ukb_snp_ind_supp,
  by = "eid",
  all = TRUE
) %>% 
  extract(, c("eid", iv_info$All_sites$affy)) %>% 
  set_colnames(c("eid", iv_info$All_sites$SNP))

snp_score <- data.frame(eid = snp_ind$eid) %>% 
  {
    for (i in iv_info$All_sites$SNP) {
      eff <- subset(
        iv_info$All_sites, 
        SNP == i, 
        "effect_allele.exposure",
        drop = TRUE
      )
      oth <- subset(
        iv_info$All_sites,
        SNP == i,
        "other_allele.exposure",
        drop = TRUE
      )
      homo_eff <- paste(eff, eff)
      homo_oth <- paste(oth, oth)
      heto <- c(paste(eff, oth), paste(oth, eff))
      
      col <- snp_ind[, i]
      for (j in seq_along(col)) {
        if (!is.na(col[j])) {
          if (col[j] == homo_oth) {
            col[j] <- 0
          } else if (col[j] == heto[1]) {
            col[j] <- 1
          } else if (col[j] == heto[2]) {
            col[j] <- 1
          } else if (col[j] == homo_eff) {
            col[j] <- 2
          } else {
            col[j] <- NA
          }
        }
      }
     .[, i] <- as.numeric(col)
    }
    .
  }

prs_all <- data.frame(
  eid = snp_score$eid,
  prs_u = rowSums(snp_score[, -1]),
  prs_w = (function() {
    if (identical(names(snp_score[, -1]), iv_info$All_sites$SNP)) {
      a <- as.matrix(snp_score[, -1])
      b <- matrix(iv_info$All_sites$beta.exposure, ncol = 1)
    } else {
      stop("Not identical in names!")
    }
    a %*% b
  })()
) %>% 
  merge(UKb_baseline[, c("eid", "platelet")], by = "eid", all.x = TRUE) %>% 
  na.omit()

prs_without_cancer <- subset(
  prs_all,
  !is.element(prs_all$eid, whole_cancer_data$All_sites$eid)
)

gg_prs(prs_all, "prs_u", "all participants")
