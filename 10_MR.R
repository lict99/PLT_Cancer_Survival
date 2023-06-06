
# env settings ------------------------------------------------------------

library(magrittr)
library(ieugwasr)
library(TwoSampleMR)
library(LDlinkR)
library(openxlsx)

load("09/Cox_regression_of_SNP.RData")
pc_snps <- read.xlsx("src/SNP/PC_SNPs.xlsx", 1)

dir.create("10", FALSE)

# data arrangement --------------------------------------------------------

# find allele frequency of EUR
if (file.exists("10/allele_frequency_of_EUR_in_LDlink.RData")) {
  load("10/allele_frequency_of_EUR_in_LDlink.RData")
} else {
  snp_eaf <- sapply(
    union(pc_snps$SNP, Cox_snp$All_sites$snp),
    function(x) {
      message("finding haplotype for ", x)
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
  save(snp_eaf, file = "10/allele_frequency_of_EUR_in_LDlink.RData")
}

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

# find proxies
snp_need_prx <- setdiff(pc_snps_ieu$SNP, Cox_snp$All_sites$snp)

if (file.exists("10/snp_proxies.RData")) {
  load("10/snp_proxies.RData")
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
  save(snp_prx, file = "10/snp_proxies.RData")
}

# select proxies
snp_prx_fil <- lapply(
  snp_prx,
  function(x) {
    subset(x, R2 >= 0.8) %>% 
      subset(is.element(RS_Number, Cox_snp$All_sites$snp)) %>% 
      extract(which.max(use_series(., R2)), )
  }
)

# substitute proxies
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

# impute EAF
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
      return(as.numeric(f))
    }
  )
)

pc_snps_exp <- format_data(
  dat = pc_snps_eaf,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "effect.allele",
  other_allele_col = "reference.allele",
  pval_col = "P",
  samplesize_col = "n",
  eaf_col = "eaf"
) %>% 
  clump_data(pop = "EUR")

css_snp_out <- lapply(
  Cox_snp,
  function(x) {
    if (any(x[, "p"] <= 5e-08, na.rm = TRUE)) {
      message("p-value <= 5E-08 in outcome!")
    }
    format_data(
      dat = x,
      type = "outcome",
      phenotype_col = "phenotype",
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
)

mr_data <- lapply(
  css_snp_out,
  function(x) {
    harmonise_data(
      exposure_dat = pc_snps_exp,
      outcome_dat = x,
      action = 2
    )
  }
)

# MR ----------------------------------------------------------------------

mr_res <- lapply(
  mr_data,
  function(x) {
    try(
      mr(x) %>% 
        generate_odds_ratios()
    )
  }
)

# data saving -------------------------------------------------------------

save(mr_data, file = "10/MR_harmonised_data.RData")
save(mr_res, file = "10/MR_results.RData")
