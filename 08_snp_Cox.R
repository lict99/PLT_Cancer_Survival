# env settings -----------------------------------------------------------------

library(magrittr)
library(survival)
library(openxlsx)
library(data.table)

load("01/whole_cancer_data_for_Cox.RData")
load("00/cancer_ICD_codes_with_attr.RData")

source("functions/Cox_regression.R", local = TRUE)

ukb_snp_ind <- fread("src/SNP/UKb_snp_individual.csv", data.table = FALSE)
ukb_snp_ind_supp <- fread(
  "src/SNP/UKb_snp_individual_supp.csv",
  data.table = FALSE
)
ukb_snp_sum <- fread("src/SNP/UKb_snp_summary.txt", data.table = FALSE)

pc_snp <- read.xlsx("src/SNP/PC_SNPs.xlsx", 1)
pc_snp_prx <- read.xlsx("src/SNP/query_snp.xlsx", 1)

dir.create("08", FALSE)

# data preprocessing -----------------------------------------------------------

## transform SNP summary information
snp_smr <- ukb_snp_sum %>%
  transform(
    affy = paste0("affy", affy_id),
    snp = paste0("rs", rs_id)
  ) %>%
  subset(subset = is.element(snp, union(pc_snp[[1]], pc_snp_prx[[1]]))) %>%
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
        switch(x[["ref"]],
          A = x[["allele_a"]],
          B = x[["allele_b"]],
          NA_character_
        )
      }
    ),
    eff_allele = apply(
      ., 1,
      function(x) {
        switch(x[["ref"]],
          A = x[["allele_b"]],
          B = x[["allele_a"]],
          NA_character_
        )
      }
    )
  ) %>%
  set_rownames(extract2(., "snp"))


## extract all patients with cancer
## and corresponding SNP information

genotype_char <- merge(ukb_snp_ind, ukb_snp_ind_supp, by = "eid", all = TRUE)

## substitute genotype with 0, 1, 2
genotype_num <- apply(
  snp_smr, 1,
  function(x) {
    homo_ref <- paste(x[["ref_allele"]], x[["ref_allele"]])
    homo_eff <- paste(x[["eff_allele"]], x[["eff_allele"]])
    hete <- c(
      paste(x[["ref_allele"]], x[["eff_allele"]]),
      paste(x[["eff_allele"]], x[["ref_allele"]])
    )
    char_col <- genotype_char[[x[["affy"]]]]
    num_col <- rep(NA_real_, times = length(char_col))
    num_col[char_col == homo_ref] <- 0
    num_col[char_col %in% hete] <- 1
    num_col[char_col == homo_eff] <- 2
    num_col
  }
) %>%
  as.data.frame() %>%
  cbind(data.frame(eid = genotype_char[["eid"]]), .)

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

cancer_eid <- Cox_data[["OS"]][["All_sites"]][["eid"]]
genotype_num_cancer <- subset(genotype_num, eid %in% cancer_eid)

## merge baseline data and SNP data
data_Cox_snp <- lapply(
  Cox_data, function(x) {
    lapply(
      x, function(y) {
        merge(y, genotype_num_cancer, by = "eid", all.x = TRUE)
      }
    )
  }
)

# Cox regression of SNP --------------------------------------------------------

## compute the association between SNPs and cancer survival
## using UKB individual data
Cox_snp <- lapply(
  data_Cox_snp, function(x) {
    lapply(
      x, function(y) {
        snp_list <- apply(
          snp_smr, 1, function(z) {
            snp <- z[["snp"]]
            if ("sex" %in% colnames(y)) {
              fml <- as.formula(
                paste("Surv(fu_time, fu_event == 1) ~ age + sex +", snp)
              )
              fit <- coxph(
                fml,
                data = y,
                singular.ok = TRUE,
                x = TRUE
              )
            } else {
              fml <- as.formula(
                paste("Surv(fu_time, fu_event == 1) ~ age +", snp)
              )
              fit <- coxph(
                fml,
                data = y,
                singular.ok = TRUE,
                x = TRUE
              )
            }
            smr <- summary(fit)
            data.frame(
              snp = snp,
              ref = z[["ref_allele"]],
              eff = z[["eff_allele"]],
              coef = smr[["coefficients"]][snp, "coef"],
              se = smr[["coefficients"]][snp, "se(coef)"],
              p = smr[["coefficients"]][snp, "Pr(>|z|)"],
              n = fit[["n"]],
              eaf = (sum(fit[["x"]][, snp])) / (fit[["n"]] * 2)
            )
          },
          simplify = FALSE
        )
        Reduce(rbind, snp_list)
      }
    )
  }
)

# data saving ------------------------------------------------------------------

save(snp_smr, file = "08/pc_snp_summary.RData", compress = FALSE)
save(genotype_char, file = "08/char_genotype.RData", compress = FALSE)
save(Cox_snp, file = "08/Cox_regression_of_SNP.RData", compress = FALSE)
