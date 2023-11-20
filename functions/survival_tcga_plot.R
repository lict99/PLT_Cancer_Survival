# KM curve plot with HR and P-value ----

#' @param gene `char` colname of a gene
#' @param data `data.frame` a data.frame of survival data
#' @param event `char` colname of survival events
#' @param fu_time `char` colname of follow-up time
#' @param eventcode `vect` code of survival events
#' @param vars `char` colnames of covariates
geom_survival <- function(
    gene,
    data,
    event,
    fu_time,
    eventcode = 1,
    vars = c("age", "gender")) {
  cox_smr <- coxph(
    as.formula(
      paste(
        paste0("Surv(", fu_time, " ,", event, " == ", eventcode, ")"),
        paste(c(gene, vars), collapse = " + "),
        sep = " ~ "
      )
    ),
    data
  ) %>% summary()
  pval <- cox_smr[["coefficients"]] %>%
    extract(grepl(gene, rownames(.)), "Pr(>|z|)")
  hr <- cox_smr[["conf.int"]] %>%
    extract(grepl(gene, rownames(.)), "exp(coef)") %>%
    sprintf("%.2f", .)
  fit <- surv_fit(
    as.formula(
      paste(
        paste0("Surv(", fu_time, " ,", event, " == ", eventcode, ")"),
        gene,
        sep = " ~ "
      )
    ),
    data = data
  )
  p1 <- ggsurvplot(
    fit,
    palette = c("#0072B5FF", "#BC3C29FF"),
    conf.int = TRUE,
    pval = FALSE,
    risk.table = TRUE,
    risk.table.y.text.col = TRUE,
    xlab = "Follow-up time (years)",
    legend = "top",
    legend.title = paste(gene, "expression"),
    legend.labs = levels(data[, gene]),
    alpha = 0.8,
    censor = FALSE
  )
  if (pval < 0.001) {
    p_ann <- "<0.001"
  } else {
    p_ann <- paste0("=='", sprintf("%.3f", pval), "'")
  }
  p1$plot <- p1$plot +
    annotate(
      "text",
      x = 0,
      y = 0.2,
      label = paste0(
        "paste(HR=='", hr, "',';',", "italic('p')", p_ann, ")"
      ),
      parse = TRUE,
      hjust = 0
    ) +
    theme(axis.title.x = element_blank())
  p1
}

# point plot of survival analysis ----

geom_point_surv <- function(data) {
  p1 <- ggplot(
    data,
    aes(
      x = HR,
      y = -log10(fdr)
    )
  ) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey") +
    geom_point(aes(color = eqtl)) +
    geom_text_repel(
      aes(
        label = ifelse(eqtl == "yes" & fdr < 0.05, gene, NA)
      ),
      seed = 1,
      nudge_y = -0.1,
      size = 3,
      na.rm = TRUE
    ) +
    annotate(
      "text",
      x = max(data$HR),
      y = -log10(0.05),
      label = "'FDR'=='0.05'",
      size = 2,
      parse = TRUE,
      hjust = 1,
      vjust = -0.5
    ) +
    facet_grid(. ~ survival) +
    scale_color_manual(
      labels = c("No", "Yes"),
      values = c("#6F99ADFF", "#BC3C29FF")
    ) +
    labs(
      color = "eQTL",
      y = expression(paste(-log[10], "FDR")),
      x = "Hazard ratio"
    ) +
    theme(
      legend.position = "top",
      panel.grid = element_blank()
    )
  p1
}
