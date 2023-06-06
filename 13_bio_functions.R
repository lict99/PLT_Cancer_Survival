# env settings

library(magrittr)
library(LDlinkR)
library(rsnps)
library(clusterProfiler)

load("10/MR_harmonised_data.RData")

dir.create("13", FALSE)

# function query

if (file.exists("13/SNP_trait.RData")) {
  load("13/SNP_trait.RData")
} else {
  trait <- LDtrait(
    snps = mr_data$All_sites$SNP,
    pop = "EUR",
    token = Sys.getenv("LDLINK_TOKEN"),
    genome_build = "grch37"
  )
  save(trait, file = "13/SNP_trait.RData")
}

if (file.exists("13/ncbi_snp.RData")) {
  load("13/ncbi_snp.RData")
} else {
  ncbi_snp <- ncbi_snp_query(mr_data$All_sites$SNP)
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
  )

go <- enrichGO(
  genelist$ENTREZID,
  "org.Hs.eg.db",
  ont = "ALL",
  readable = TRUE
)

kegg <- enrichKEGG(
  genelist$ENTREZID,
  "hsa"
) %>% 
  setReadable(
    "org.Hs.eg.db",
    keyType = "ENTREZID"
  )

mkegg <- enrichMKEGG(
  genelist$ENTREZID
)

enrichWP(
  genelist$ENTREZID,
  "Homo sapiens"
)

dotplot(go)
dotplot(kegg)
barplot(go)
barplot(kegg)
