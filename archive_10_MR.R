
# env settings ------------------------------------------------------------

library(magrittr)
library(TwoSampleMR)
library(openxlsx)
library(rsnps)
library(LDlinkR)

load("09/Cox_regression_of_SNP.RData")
pc_snps <- read.xlsx("src/SNP/PC_SNPs.xlsx", 1)

dir.create("10", FALSE)

# data arrangement --------------------------------------------------------

if (grepl("ncbi_pc", list.files("10", "RData"))) {
  load("10/ncbi_pc_2023-05-03.RData")
} else {
  ncbi_pc <- ncbi_snp_query(snp_pc$SNP)
  save(ncbi_pc, file = paste0("10/ncbi_pc_", Sys.Date(), ".RData"))
}

if (all(ncbi_pc$rsid == pc_snps$SNP)) {
  snp_pc <- transform(
    pc_snps,
    chr = ncbi_pc$chromosome,
    pos = ncbi_pc$bp,
    phenotype = "PC",
    eaf = (function() {
      res <- numeric()
      for (i in 1:nrow(ncbi_pc)) {
        mmaf <- ncbi_pc[[16]][[i]] %>% 
          subset(study == "dbGaP_PopFreq") %>% 
          subset(MAF != 0)
        eff <- pc_snps[i, 2]
        if (mmaf[, "Minor"] == eff) {
          maf <- mmaf[, "MAF"]
        } else if (mmaf[, "ref_seq"] == eff) {
          maf <- 1 - (mmaf[, "MAF"])
        } else {
          maf <- NA_real_
        }
        res <- c(res, maf)
      }
      return(res)
    })()
  ) %>% 
    format_data(
      type = "exposure",
      phenotype_col = "phenotype",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "SE",
      eaf_col = "eaf",
      effect_allele_col = "effect.allele",
      other_allele_col = "reference.allele",
      pval_col = "P",
      samplesize_col = "n",
      chr_col = "chr",
      pos_col = "pos"
    )
}

snp_surv <- lapply(
  Cox_snp,
  function(x) {
    transform(
      x,
      phenotype = "CSS"
    ) %>% 
      format_data(
        type = "outcome",
        phenotype_col = "phenotype",
        snp_col = "snp",
        beta_col = "coef",
        se_col = "se",
        eaf_col = "maf",
        effect_allele_col = "eff",
        other_allele_col = "ref",
        pval_col = "p"
      )
  }
)

##
##
##
palindromic_snp <- harmonise_data(snp_pc, snp_surv$All_sites, action = 2) %>% 
  extract(
    use_series(., palindromic) == TRUE & use_series(., ambiguous) == TRUE, 
    "SNP"
  )

absent_snp <- setdiff(snp_pc$SNP, snp_surv$All_sites$SNP)

proxy <- list()
for (i in c(palindromic_snp, absent_snp)) {
  message(paste("Finding proxies for", i))
  tryCatch(
    {
      proxy[[i]] <- LDproxy(
        snp = i,
        pop = "EUR",
        token = Sys.getenv("LDLINK_TOKEN"),
        genome_build = "grch37"
      )
    },
    error = function(e) {
      NA
    }
  )
}
proxy_fil <- lapply(proxy, function(x) subset(x, R2 >= 0.8) %>% try(silent = TRUE))
write.xlsx(c(list(palindromic_snp = palindromic_snp, absent_snp = absent_snp), proxy_fil),
           file = "src/SNP/proxies_snp.xlsx", overwrite = FALSE)
a <- read.xlsx("src/SNP/query_snp.xlsx", 1, colNames = FALSE)
b <- setdiff(a$X1, snp_surv$All_sites$SNP)
write.xlsx(data.frame(snp = b), file = "src/SNP/query_snp_fil.xlsx", overwrite = FALSE)
##
##

mr_data <- harmonise_data(snp_pc, snp_surv$All_sites, action = 2)


mr_res <- mr(mr_data) %>% generate_odds_ratios()
mr_heto <- mr_heterogeneity(mr_data)
mr_sin <- mr_singlesnp(mr_data)
mr_forest_plot(mr_sin)
mr_density_plot(mr_sin, mr_res)
mr_loo <- mr_leaveoneout(mr_data)
mr_leaveoneout_plot(mr_loo)
mr_ple <- mr_pleiotropy_test(mr_data)
mr_funnel_plot(mr_sin)
