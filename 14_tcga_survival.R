# env settings -----------------------------------------------------------------

library(magrittr)
library(openxlsx)
library(survival)
library(survminer)
library(ggplot2)
library(ggrepel)

source("functions/survival_tcga_plot.R", local = TRUE)

load("src/tcga/gene_expression.RData")
load("13/genes_by_loci.RData")

clinical <- read.xlsx("src/tcga/clinical_data.xlsx", 1)
eqtl_index <- read.xlsx("13/SNP_info_from_NCBI.xlsx")

dir.create("14", FALSE)

# combination of gene expression and clinical information ----------------------

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
) %>%
  (function(x) {
    for (i in seq_len(nrow(x))) {
      if (anyNA(c(x[i, "OS_time"], x[i, "OS"]))) {
        x[i, "OS_time"] <- x[i, "OS"] <- NA
      } else if (x[i, "OS_time"] > 10) {
        x[i, "OS_time"] <- 10
        if (x[i, "OS"] == 1) {
          x[i, "OS"] <- 0
        }
      }
    }
    x
  })() %>%
  (function(x) {
    for (i in seq_len(nrow(x))) {
      if (anyNA(c(x[i, "CSS_time"], x[i, "CSS"]))) {
        x[i, "CSS_time"] <- x[i, "CSS"] <- NA
      } else if (x[i, "CSS_time"] > 10) {
        x[i, "CSS_time"] <- 10
        if (x[i, "CSS"] == 1) {
          x[i, "CSS"] <- 0
        }
      }
    }
    x
  })()

pan_surv_df <- pan_data %>%
  (function(x) {
    df <- data.frame()
    for (g in genes) {
      os_fit <- coxph(
        as.formula(
          paste(
            "Surv(OS_time, OS)",
            paste(c(g, "age", "gender"), collapse = " + "),
            sep = " ~ "
          )
        ),
        data = x
      ) %>% summary()
      css_fit <- coxph(
        as.formula(
          paste(
            "Surv(CSS_time, CSS)",
            paste(c(g, "age", "gender"), collapse = " + "),
            sep = " ~ "
          )
        ),
        data = x
      ) %>% summary()
      hr_os <- os_fit[["conf.int"]] %>%
        extract(grepl(g, rownames(.)), "exp(coef)") %>%
        sprintf("%.2f", .)
      hr_os_l95 <- os_fit[["conf.int"]] %>%
        extract(grepl(g, rownames(.)), "lower .95") %>%
        sprintf("%.2f", .)
      hr_os_u95 <- os_fit[["conf.int"]] %>%
        extract(grepl(g, rownames(.)), "upper .95") %>%
        sprintf("%.2f", .)
      p_os <- os_fit[["coefficients"]] %>%
        extract(grepl(g, rownames(.)), "Pr(>|z|)")
      hr_css <- css_fit[["conf.int"]] %>%
        extract(grepl(g, rownames(.)), "exp(coef)") %>%
        sprintf("%.2f", .)
      hr_css_l95 <- css_fit[["conf.int"]] %>%
        extract(grepl(g, rownames(.)), "lower .95") %>%
        sprintf("%.2f", .)
      hr_css_u95 <- css_fit[["conf.int"]] %>%
        extract(grepl(g, rownames(.)), "upper .95") %>%
        sprintf("%.2f", .)
      p_css <- css_fit[["coefficients"]] %>%
        extract(grepl(g, rownames(.)), "Pr(>|z|)")
      df_g <- rbind(
        data.frame(
          gene = g,
          survival = "OS",
          HR = hr_os,
          lower.95 = hr_os_l95,
          upper.95 = hr_os_u95,
          pval = p_os
        ),
        data.frame(
          gene = g,
          survival = "CSS",
          HR = hr_css,
          lower.95 = hr_css_l95,
          upper.95 = hr_css_u95,
          pval = p_css
        )
      )
      df <- rbind(df, df_g)
    }
    df[, "fdr"] <- p.adjust(df[, "pval"], method = "fdr")
    df[, "eqtl"] <- ifelse(
      is.element(df[, "gene"], eqtl_index[eqtl_index$eqtl == "yes", "gene"]),
      "yes",
      "no"
    )
    df
  })()

# plotting ---------------------------------------------------------------------

## scatter plot of survival

surv_plot <- geom_point_surv(
  data = transform(
    pan_surv_df,
    HR = as.numeric(HR),
    eqtl = ifelse(eqtl == "yes" & fdr < 0.05, "yes", "no"),
    survival = factor(
      survival,
      levels = c("OS", "CSS"),
      labels = c("Overall survival", "Cancer-specific survival")
    )
  )
)

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

# results saving ---------------------------------------------------------------

write.xlsx(
  pan_surv_df,
  "14/pancancer_survival_summary.xlsx",
  TRUE
)

# plots saving -----------------------------------------------------------------

pdf(paste0("14/OS.pdf"), 7, 7)
print(pan_survplot_os)
dev.off()

pdf(paste0("14/CSS.pdf"), 7, 7)
print(pan_survplot_css)
dev.off()

ggsave(
  surv_plot,
  filename = "14/survival_scatter.pdf",
  width = 7,
  height = 5
)

# data saving ------------------------------------------------------------------

tpm4_os <- pan_survplot_os[["TPM4"]]
save(tpm4_os, file = "14/km_tpm4_os.RData", compress = FALSE)
