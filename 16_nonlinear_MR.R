# env settings ----

library(magrittr)
library(nlmr)
library(ggplot2)
library(metafor)
library(patchwork)

load("11/PRS_by_cancer.RData")
load("00/cancer_names.RData")

source("functions/nonlinear_MR_plot.R", local = TRUE)

dir.create("16", FALSE)

# ----

nlmr_plot <- lapply(
  prs_cancer,
  function(x) {
    lapply(
      x[c(1)],
      function(y) {
        for (i in seq_len(nrow(y))) {
          if (y[i, "fu_time"] > 365.25 * 5) {
            y[i, "fu_time"] <- 365.25 * 5
            if (y[i, "fu_event"] == 1) {
              y[i, "fu_event"] <- 0
            }
          }
        }
        set.seed(1)
        res_w <- fracpoly_mr2(
          y = y[, "fu_event"],
          x = y[, "platelet"],
          g = y[, "prs_w"],
          d = "both",
          covar = if ("sex" %in% colnames(y)) {
            y[, c("age", "sex")]
          } else {
            y[, "age", drop = FALSE]
          },
          family = "binomial",
          xpos = 0.5,
          method = "REML",
          ref = 400,
          fig = TRUE
        )
        set.seed(1)
        res_u <- fracpoly_mr2(
          y = y[, "fu_event"],
          x = y[, "platelet"],
          g = y[, "prs_u"],
          d = "both",
          covar = if ("sex" %in% colnames(y)) {
            y[, c("age", "sex")]
          } else {
            y[, "age", drop = FALSE]
          },
          family = "binomial",
          xpos = 0.5,
          method = "REML",
          ref = 400,
          fig = TRUE
        )
        list(
          prs_w = res_w$figure +
            labs(
              subtitle = "Using weighted PRS as the instrumental variable"
            ),
          prs_u = res_u$figure +
            labs(
              subtitle = "Using unweighted PRS as the instrumental variable"
            )
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
          function(y) y + labs(title = paste("Five-year", ii, "in", jj))
        )
      }
    }
    x
  })() %>%
  unlist(FALSE) %>%
  unlist(FALSE) %>%
  {
    wrap_plots(., guides = "collect") +
      plot_annotation(tag_levels = "A") &
      theme(legend.position = "bottom")
  }

# plots saving ----

ggsave(
  "16/nonlinear_MR_plot.pdf",
  nlmr_plot,
  height = 10,
  width = 10,
  device = "pdf"
)
