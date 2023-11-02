# env settings ----

library(magrittr)
library(nlmr)
library(ggplot2)
library(metafor)

load("11/PRS_by_cancer.RData")
load("00/cancer_names.RData")

source("functions/nonlinear_MR_plot.R")

dir.create("16", FALSE)

# ----

nlmr_plot <- lapply(
  prs_cancer,
  function(x) {
    lapply(
      x[c(1, 2)],
      function(y) {
        for (i in seq_len(nrow(y))) {
          if (y[i, "fu_time"] > 365.25 * 5) {
            y[i, "fu_time"] <- 365.25 * 5
            if (y[i, "fu_event"] == 1) {
              y[i, "fu_event"] <- 0
            }
          }
        }
        res_w <- fracpoly_mr2(
          y = y[, "fu_event"],
          x = y[, "platelet"],
          g = y[, "prs_w"],
          covar = if ("sex" %in% colnames(y)) {
            y[, c("age", "sex")]
          } else {
            y[, "age", drop = FALSE]
          },
          family = "binomial",
          xpos = 0.5,
          method = "REML",
          ref = median(y[, "platelet"], na.rm = TRUE),
          fig = TRUE
        )
        res_u <- fracpoly_mr2(
          y = y[, "fu_event"],
          x = y[, "platelet"],
          g = y[, "prs_u"],
          covar = if ("sex" %in% colnames(y)) {
            y[, c("age", "sex")]
          } else {
            y[, "age", drop = FALSE]
          },
          family = "binomial",
          xpos = 0.5,
          method = "REML",
          ref = median(y[, "platelet"], na.rm = TRUE),
          fig = TRUE
        )
        list(
          prs_w = res_w$figure + labs(subtitle = "Weighted PRS"),
          prs_u = res_u$figure + labs(subtitle = "Unweighted PRS")
        )
      }
    )
  }
) %>%
  (function(x) {
    for (i in names(x)) {
      for (j in names(x[[i]])) {
        ii <- switch(i,
          OS = "overall survival",
          CSS = "cancer-specific survival",
        )
        jj <- cancer_names[[j]] %>% tolower()
        x[[i]][[j]] <- lapply(
          x[[i]][[j]],
          function(y) y + labs(title = paste("5-year", ii, "in", jj))
        )
      }
    }
    x
  })()

# plots saving ----

pdf("16/nonlinear_MR_plot.pdf", width = 8, height = 6)
print(nlmr_plot)
dev.off()
