# env settings -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)

load("01/whole_cancer_data_for_Cox.RData")
load("08/char_genotype.RData")
load("08/pc_snp_summary.RData")
load("09/MR_harmonised_data.RData")
load("00/cancer_ICD_codes_with_attr.RData")
load("00/cancer_names.RData")

source("functions/Cox_regression.R", local = TRUE)
source("functions/PRS_plot.R", local = TRUE)

dir.create("11", FALSE)

# PRS calculation --------------------------------------------------------------

## extract IV information
iv_info <- mr_data[["OS"]][["All_sites"]] %>%
  subset(mr_keep == TRUE) %>%
  merge(
    snp_smr[, c("affy", "snp")],
    by.x = "SNP",
    by.y = "snp",
    all.x = TRUE
  ) %>%
  set_rownames(extract2(., "SNP"))

## calculate SNP effect allele score
genotype_num <- apply(
  iv_info, 1,
  function(x) {
    ref_allele <- x[["other_allele.exposure"]]
    eff_allele <- x[["effect_allele.exposure"]]
    homo_ref <- paste(ref_allele, ref_allele)
    homo_eff <- paste(eff_allele, eff_allele)
    hete <- c(
      paste(ref_allele, eff_allele),
      paste(eff_allele, ref_allele)
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

## unweighted and weighted PRS for various cancers
prs <- data.frame(
  eid = genotype_num[["eid"]],
  prs_u = rowSums(genotype_num[, -1]),
  prs_w = (function() {
    stopifnot(identical(names(genotype_num[, -1]), iv_info[["SNP"]]))
    a <- as.matrix(genotype_num[, -1])
    b <- matrix(iv_info[["beta.exposure"]], ncol = 1)
    a %*% b
  })()
)

## PRS of cancer patients
prs_cancer <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = c("eid", "platelet", "fu_time", "fu_event", "age", "sex"),
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

# plotting ---------------------------------------------------------------------

## weighted PRS plots for different cancer types
prs_w_plot <- mapply(
  function(x, y) {
    geom_prs(
      data = x,
      annotate = cancer_names[[y]],
      x_col = "prs_w"
    )
  },
  prs_cancer[["OS"]],
  names(prs_cancer[["OS"]]),
  SIMPLIFY = FALSE
)


## unweighted PRS plots for different cancer types
prs_u_plot <- mapply(
  function(x, y) {
    geom_prs(
      data = x,
      annotate = cancer_names[[y]],
      x_col = "prs_u"
    )
  },
  prs_cancer[["OS"]],
  names(prs_cancer[["OS"]]),
  SIMPLIFY = FALSE
)

# plots saving -----------------------------------------------------------------

pdf("11/weighted_PRS_in_cancer.pdf", width = 9.5, height = 9.5)
print(prs_w_plot)
dev.off()

pdf("11/unweighted_PRS_in_cancer.pdf", width = 9.5, height = 9.5)
print(prs_u_plot)
dev.off()

# data saving ------------------------------------------------------------------

save(iv_info, file = "11/IV_info.RData", compress = FALSE)
save(prs_cancer, file = "11/PRS_by_cancer.RData", compress = FALSE)
save(prs, file = "11/PRS.RData", compress = FALSE)
