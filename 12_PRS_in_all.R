# env settings ----

library(magrittr)
library(survival)
library(ggplot2)
library(patchwork)
library(data.table)

load("11/IV_info.RData")
load("01/UKB_all_info.RData")

source("functions/PRS_plot.R")

ukb_snp_ind <- fread("src/SNP/UKb_snp_individual.csv", data.table = FALSE)
ukb_snp_ind_supp <- fread(
  "src/SNP/UKb_snp_individual_supp.csv",
  data.table = FALSE
)

dir.create("12", FALSE)

# calculation ----

## extract individual data of valid genetic instruments
snp_ind <- merge(
  ukb_snp_ind,
  ukb_snp_ind_supp,
  by = "eid",
  all = TRUE
) %>%
  extract(, c("eid", iv_info$affy)) %>%
  set_colnames(c("eid", iv_info$SNP))

## calculate SNP effect allele score
snp_score <- data.frame(eid = snp_ind$eid) %>%
  (function(score) {
    for (i in iv_info$SNP) {
      eff <- subset(
        iv_info,
        SNP == i,
        "effect_allele.exposure",
        drop = TRUE
      )
      oth <- subset(
        iv_info,
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
      score[, i] <- as.numeric(col)
    }
    score
  })()

## PRS for all participants in the UK biobank
prs_all <- data.frame(
  eid = snp_score$eid,
  prs_u = rowSums(snp_score[, -1]),
  prs_w = (function() {
    if (identical(names(snp_score[, -1]), iv_info$SNP)) {
      a <- as.matrix(snp_score[, -1])
      b <- matrix(iv_info$beta.exposure, ncol = 1)
    } else {
      stop("Not identical in names!")
    }
    a %*% b
  })()
) %>%
  merge(
    UKb_baseline[, c("eid", "platelet")],
    by = "eid",
    all.x = TRUE
  ) %>%
  na.omit()

# plotting ----

## weighted PRS plot
prs_w_plot <- geom_prs(
  data = prs_all,
  annotate = "All participants in UK Biobank",
  x_col = "prs_w"
)

## unweighted PRS plot
prs_u_plot <- geom_prs(
  data = prs_all,
  annotate = "All participants in UK Biobank",
  x_col = "prs_u"
)

# plots saving ----

ggsave(
  "12/weighted_PRS_in_all.pdf",
  prs_w_plot,
  width = 13.4,
  height = 9.91,
  units = "in",
  device = grDevices::cairo_pdf
)

ggsave(
  "12/unweighted_PRS_in_all.pdf",
  prs_u_plot,
  width = 13.4,
  height = 9.91,
  units = "in",
  device = grDevices::cairo_pdf
)
