# env settings -----------------------------------------------------------------

library(magrittr)
library(ieugwasr)
library(TwoSampleMR)
library(LDlinkR)
library(openxlsx)

load("00/cancer_names.RData")
load("08/Cox_regression_of_SNP.RData")

pc_snps <- read.xlsx("src/SNP/PC_SNPs.xlsx", 1)
pheno_df <- read.xlsx("src/SNP/phenoscanner.xlsx", 1)
confounder <- read.xlsx("src/SNP/phenoscanner.xlsx", 2)

dir.create("09", FALSE)

# data arrangement -------------------------------------------------------------

## SNPs in Cox regression
snps_cox <- Cox_snp[["OS"]][["All_sites"]][["snp"]]

## find allele frequencies of EUR using LDlink
if (file.exists("src/cache/09_allele_frequency_of_EUR_in_LDlink.RData")) {
  load("src/cache/09_allele_frequency_of_EUR_in_LDlink.RData")
} else {
  snp_eaf <- sapply(
    union(pc_snps$SNP, snps_cox),
    function(x) {
      message("Finding haplotype for ", x)
      try(
        LDhap(
          snps = x,
          pop = "EUR",
          token = Sys.getenv("LDLINK_TOKEN"),
          table_type = "haplotype",
          genome_build = "grch37"
        )
      )
    },
    simplify = FALSE
  )
  save(snp_eaf, file = "src/cache/09_allele_frequency_of_EUR_in_LDlink.RData")
}

## linkage disequilibrium test
if (file.exists("src/cache/09_pc_snp_ieu.RData")) {
  load("src/cache/09_pc_snp_ieu.RData")
} else {
  pc_snps_ieu <- subset(
    pc_snps,
    is.element(SNP, ld_reflookup(rsid = SNP, pop = "EUR"))
  ) %>%
    format_data(
      type = "exposure",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "SE",
      effect_allele_col = "effect.allele",
      other_allele_col = "reference.allele",
      pval_col = "P",
      samplesize_col = "n"
    ) %>%
    clump_data(pop = "EUR")
  save(pc_snps_ieu, file = "src/cache/09_pc_snp_ieu.RData")
}

## find proxies
snp_need_prx <- setdiff(
  pc_snps_ieu$SNP,
  snps_cox
)

if (file.exists("src/cache/09_snp_proxies.RData")) {
  load("src/cache/09_snp_proxies.RData")
} else {
  snp_prx <- sapply(
    snp_need_prx,
    function(x) {
      message("finding proxies for ", x)
      LDproxy(
        snp = x,
        pop = "EUR",
        r2d = "r2",
        token = Sys.getenv("LDLINK_TOKEN"),
        genome_build = "grch37"
      )
    },
    simplify = FALSE
  )
  save(snp_prx, file = "src/cache/09_snp_proxies.RData")
}

## select proxies
snp_prx_fil <- lapply(
  snp_prx,
  function(x) {
    subset(x, R2 >= 0.8) %>%
      subset(is.element(RS_Number, snps_cox)) %>%
      extract(which.max(use_series(., R2)), )
  }
)

## substitute proxies
pc_snps_prx <- pc_snps %>%
  subset(is.element(SNP, pc_snps_ieu$SNP))
for (i in snp_need_prx) {
  prx_sub <- snp_prx_fil[[i]][, "RS_Number"]
  a_info <- snp_prx_fil[[i]][, "Correlated_Alleles"] %>%
    strsplit("", fixed = TRUE) %>%
    extract2(1)
  if (length(a_info) == 7) {
    if (all(pc_snps[pc_snps$SNP == i, c(2, 3)] == a_info[c(1, 5)])) {
      prx_sub[c(2, 3)] <- a_info[c(3, 7)]
    } else if (all(pc_snps[pc_snps$SNP == i, c(2, 3)] == a_info[c(5, 1)])) {
      prx_sub[c(2, 3)] <- a_info[c(7, 3)]
    } else {
      stop("Mismatch allele!")
    }
  } else {
    stop("Not SNP proxy!")
  }
  pc_snps_prx[pc_snps_prx$SNP == i, ][1:3] <- prx_sub
}

## impute effect allele frequency (EAF)
pc_snps_eaf <- transform(
  pc_snps_prx,
  eaf = apply(
    pc_snps_prx, 1,
    function(x) {
      eaf_df <- snp_eaf[[x[1]]]
      if (x[2] == eaf_df[1, 1]) {
        f <- eaf_df[1, 3]
      } else if (x[2] == eaf_df[2, 1]) {
        f <- eaf_df[2, 3]
      } else {
        stop("Mismatch allele!")
      }
      as.numeric(f)
    }
  )
)

## second round of LD clumping
if (file.exists("src/cache/09_pc_snp_exp_all.RData")) {
  load("src/cache/09_pc_snp_exp_all.RData")
} else {
  pc_snps_exp_all <- format_data(
    dat = transform(
      pc_snps_eaf,
      Phenotype = "platelet counts",
      id = "PLT10^9/L"
    ),
    type = "exposure",
    snp_col = "SNP",
    phenotype_col = "Phenotype",
    beta_col = "Effect",
    se_col = "SE",
    effect_allele_col = "effect.allele",
    other_allele_col = "reference.allele",
    pval_col = "P",
    samplesize_col = "n",
    eaf_col = "eaf",
    id_col = "id"
  ) %>%
    clump_data(pop = "EUR")
  save(pc_snps_exp_all, file = "src/cache/09_pc_snp_exp_all.RData")
}

## confounding SNPs
cfd_trait <- subset(confounder, mark == 1)$trait %>% unique()
cfd_snp <- subset(pheno_df, is.element(trait, cfd_trait))$snp %>% unique()

## remove confounding SNPs
pc_snps_exp <- subset(pc_snps_exp_all, !is.element(SNP, cfd_snp))

## format data regarding the effect of SNPs on outcomes
snp_out <- list()
for (i in names(Cox_snp)) {
  for (j in names(Cox_snp[[i]])) {
    if (any(Cox_snp[[i]][[j]][, "p"] <= 5e-08, na.rm = TRUE)) {
      stop("p-value <= 5E-08 in outcome!")
    }
    snp_out[[i]][[j]] <- format_data(
      dat = transform(
        Cox_snp[[i]][[j]],
        Phenotype = i,
        id = cancer_names[[j]]
      ),
      type = "outcome",
      phenotype_col = "Phenotype",
      snp_col = "snp",
      beta_col = "coef",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "eff",
      other_allele_col = "ref",
      pval_col = "p",
      samplesize_col = "n"
    )
  }
}

## harmonize data
mr_data <- lapply(
  snp_out,
  function(x) {
    lapply(
      x,
      function(y) {
        harmonise_data(
          exposure_dat = pc_snps_exp,
          outcome_dat = y,
          action = 2
        )
      }
    )
  }
)

##
mr_data_outcome <- lapply(
  mr_data,
  function(x) {
    lapply(
      x,
      function(y) {
        df <- subset(y, mr_keep) %>%
          extract(c(
            "SNP",
            "effect_allele.outcome",
            "other_allele.outcome",
            "id.outcome",
            "beta.outcome",
            "se.outcome",
            "pval.outcome",
            "eaf.outcome"
          )) %>%
          transform(
            beta.outcome = sprintf("%.3f", beta.outcome),
            se.outcome = sprintf("%.3f", se.outcome),
            pval.outcome = sprintf("%.3f", pval.outcome),
            eaf.outcome = sprintf("%.3f", eaf.outcome)
          )
      }
    )
  }
)

# results saving ---------------------------------------------------------------

write.xlsx(mr_data_outcome[["OS"]], "09/MR_outcome_data_OS.xlsx", TRUE)
write.xlsx(mr_data_outcome[["CSS"]], "09/MR_outcome_data_CSS.xlsx", TRUE)

# data saving ------------------------------------------------------------------

save(mr_data, file = "09/MR_harmonised_data.RData", compress = FALSE)
save(snp_out, file = "09/snp_outcomes.RData", compress = FALSE)
