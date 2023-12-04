# env settings -----------------------------------------------------------------

library(magrittr)
library(TwoSampleMR)
library(data.table)
library(openxlsx)

load("13/eQTL.RData")
load("09/snp_outcomes.RData")

gtex <- fread(
  "src/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",
  data.table = FALSE
)

tcga_surv <- read.xlsx("14/pancancer_survival_summary.xlsx", 1)
ncbi_snp <- read.xlsx("13/SNP_info_from_NCBI.xlsx", 1)

dir.create("15", FALSE)

# calculation ------------------------------------------------------------------

eqtl_filter <- subset(
  eqtl,
  sentinel %in% ncbi_snp$rsid &
    !is.na(beta) &
    gene %in% tcga_surv[tcga_surv$fdr <= 0.05, "gene"]
)

gtex_filter <- subset(
  gtex,
  substr(gene_id, 1, 15) %in% eqtl_filter$ensgid
) %>%
  subset(
    sapply(strsplit(variant_id, "_"), function(x) x[2]) %in%
      eqtl_filter$var_pos_hg38
  ) %>%
  transform(
    ref = sapply(strsplit(variant_id, "_"), function(x) x[3]),
    alt = sapply(strsplit(variant_id, "_"), function(x) x[4]),
    gene_id = substr(gene_id, 1, 15)
  ) %>%
  merge(
    eqtl_filter[
      ,
      c("sentinel", "chr", "var_pos_hg38", "gene", "ensgid", "type")
    ],
    by.x = "gene_id",
    by.y = "ensgid"
  )

snp_expo <- format_data(
  dat = gtex_filter,
  type = "exposure",
  phenotype_col = "type",
  snp_col = "sentinel",
  beta_col = "slope",
  se_col = "slope_se",
  eaf_col = "maf", # maybe
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval_nominal",
  id_col = "gene_id",
  gene_col = "gene"
)

mr_data_eqtl <- lapply(
  snp_out,
  function(x) {
    lapply(
      x,
      function(y) {
        harmonise_data(
          exposure_dat = snp_expo,
          outcome_dat = y,
          action = 2
        )
      }
    )
  }
)

mr_res_eqtl <- lapply(
  mr_data_eqtl,
  function(x) {
    lapply(
      x,
      function(y) {
        mr(
          y,
          method_list = c(
            "mr_wald_ratio"
          )
        ) %>%
          generate_odds_ratios() %>%
          transform(
            OR_f = paste0(
              sprintf("%.3f", or),
              " (",
              sprintf("%.3f", or_lci95), "-", sprintf("%.3f", or_uci95),
              ")"
            ),
            p_f = ifelse(
              pval < 0.001,
              sprintf("%.3e", pval),
              sprintf("%.3f", pval)
            )
          )
      }
    )
  }
)

# results saving ---------------------------------------------------------------

write.xlsx(
  lapply(
    mr_res_eqtl,
    function(x) {
      Reduce(rbind, x)
    }
  ),
  "15/MR_results_eQTL.xlsx",
  TRUE
)
