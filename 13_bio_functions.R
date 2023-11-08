# env settings ----

library(magrittr)
library(rsnps)
library(clusterProfiler)
library(Qtlizer)
library(openxlsx)

load("11/IV_info.RData")

dir.create("13", FALSE)

# function query ----

if (file.exists("src/cache/13_ncbi_snp.RData")) {
  load("src/cache/13_ncbi_snp.RData")
} else {
  ncbi_snp <- ncbi_snp_query(iv_info$SNP)
  save(ncbi_snp, file = "src/cache/13_ncbi_snp.RData")
}

genelist <- ncbi_snp$gene %>%
  strsplit("/", fixed = TRUE) %>%
  unlist() %>%
  bitr(
    "SYMBOL",
    "ENTREZID",
    "org.Hs.eg.db",
    TRUE
  ) %>%
  {
    alias <- bitr(
      extract2(., "SYMBOL"),
      "SYMBOL",
      "ALIAS",
      "org.Hs.eg.db",
      FALSE
    )
    merge(., alias, by = "SYMBOL", all = TRUE)
  }

if (!file.exists("src/cache/13_QTL.RData")) {
  qtl <- get_qtls(genelist$ALIAS, corr = NA, ref_version = "hg19")
  save(qtl, file = "src/cache/13_QTL.RData")
} else {
  load("src/cache/13_QTL.RData")
}

eqtl <- subset(
  qtl,
  p <= 5e-08 &
    type == "eQTL" &
    grepl("blood", tissue, ignore.case = TRUE) &
    !grepl("cell", tissue, ignore.case = TRUE) &
    is.element(gene, genelist$ALIAS) &
    source == "GTEx v8"
)

ncbi_snp[, "eqtl"] <- ifelse(ncbi_snp$rsid %in% eqtl$sentinel, "yes", "no")

# data saving ----

save(genelist, file = "13/genes_by_loci.RData")
save(eqtl, file = "13/eQTL.RData")
save(ncbi_snp, file = "13/ncbi_snp.RData")
write.xlsx(ncbi_snp, "13/SNP_info_from_NCBI.xlsx", TRUE)
