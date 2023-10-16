# env settings ----

library(magrittr)
library(openxlsx)
library(survival)
library(survminer)
library(ggplot2)

source("functions/survival_tcga_plot.R")

load("src/tcga/gene_expression.RData")
load("13/genes_by_loci.RData")

clinical <- read.xlsx("src/tcga/clinical_data.xlsx", 1)

dir.create("14", FALSE)

# combination of gene expression and clinical information ----

alias_index <- tapply(genelist[, 3], genelist[, 1], c)

gene_tpm <- lapply(
  gene_expr,
  function(x, index) {
    info <- x[["geneinfo"]]
    symbol2id <- info[which(is.element(info[, "gene_name"], unlist(index))), ]
    tpm <- x[["tpm"]]
    df <- data.frame(sample = substr(colnames(tpm), 1, 12))
    for (i in names(index)) {
      tpm_i <- tpm[
        symbol2id[is.element(symbol2id[, "gene_name"], index[[i]]), "gene_id"],
      ]
      if (nrow(tpm_i) == 0L) {
        next
      } else if (nrow(tpm_i) > 1L) {
        tpm_i <- colSums(tpm_i, na.rm = TRUE)
      }
      tpm_i <- as.numeric(tpm_i)
      df[, i] <- factor(
        ifelse(tpm_i > median(tpm_i, na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High")
      )
    }
    df
  },
  index = alias_index
) %>%
  {
    pan_tpm <- data.frame()
    for (i in names(.)) {
      pan_tpm <- rbind(pan_tpm, extract2(., i))
    }
    inset2(., "pan", pan_tpm)
  }

genes <- names(gene_tpm[["pan"]]) %>%
  extract(not(equals(., "sample")))

c_clinical <- with(
  clinical,
  data.frame(
    sample = bcr_patient_barcode,
    type = type,
    age = age_at_initial_pathologic_diagnosis,
    gender = factor(gender, levels = c("MALE", "FEMALE")),
    OS = OS,
    OS_time = OS.time / 365.25,
    CSS = DSS,
    CSS_time = DSS.time / 365.25
  )
)

pan_data <- merge(
  gene_tpm[["pan"]],
  c_clinical,
  all.x = TRUE,
  by = "sample",
  sort = FALSE
)

# plotting ----

## survival plots of OS
pan_survplot_os <- sapply(
  genes,
  geom_survival,
  data = pan_data,
  event = "OS",
  fu_time = "OS_time",
  simplify = FALSE
)

## survival plots of CSS
pan_survplot_css <- sapply(
  genes,
  geom_survival,
  data = pan_data,
  event = "CSS",
  fu_time = "CSS_time",
  simplify = FALSE
)

# plot saving ----

pdf(paste0("14/OS.pdf"), 7, 7)
print(pan_survplot_os)
dev.off()

pdf(paste0("14/CSS.pdf"), 7, 7)
print(pan_survplot_css)
dev.off()
