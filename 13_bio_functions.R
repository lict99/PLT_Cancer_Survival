# env settings ----

library(magrittr)
library(rsnps)
library(clusterProfiler)
library(openxlsx)

load("11/IV_info.RData")

dir.create("13", FALSE)

# function query ----

if (file.exists("13/ncbi_snp.RData")) {
  load("13/ncbi_snp.RData")
} else {
  ncbi_snp <- ncbi_snp_query(iv_info$SNP)
  save(ncbi_snp, file = "13/ncbi_snp.RData")
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
    ensembl <- bitr(
      extract2(., "SYMBOL"),
      "SYMBOL",
      "ALIAS",
      "org.Hs.eg.db",
      FALSE
    )
    merge(., ensembl, by = "SYMBOL", all = TRUE)
  }

# go <- enrichGO(
#   genelist$ENTREZID,
#   "org.Hs.eg.db",
#   ont = "ALL",
#   readable = TRUE,
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# )

# kegg <- enrichKEGG(
#   genelist$ENTREZID,
#   "hsa",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# ) %>%
#   setReadable(
#     "org.Hs.eg.db",
#     keyType = "ENTREZID"
#   )

# mkegg <- enrichMKEGG(
#   genelist$ENTREZID,
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# )

# wp <- enrichWP(
#   genelist$ENTREZID,
#   "Homo sapiens",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# )

# # plot saving ----

# pdf("13/enrichment_dotplot.pdf")
# dotplot(go)
# dev.off()

# data saving ----

save(genelist, file = "13/genes_by_loci.RData")
write.xlsx(ncbi_snp, "13/SNP_info_from_NCBI.xlsx", TRUE)
