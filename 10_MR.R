# env settings ----

library(magrittr)
library(TwoSampleMR)
library(parallel)
library(openxlsx)
library(ggplot2)
library(patchwork)

load("00/cancer_names.RData")
load("09/MR_harmonised_data.RData")

source("functions/leave_one_out_plot.R")

dir.create("10", FALSE)

set.seed(1)

# MR ----

## MR analysis using IVW, MR-Egger, and weighted median methods
mr_res <- lapply(
  mr_data,
  function(x) {
    lapply(
      x,
      function(y) {
        mr(
          y,
          method_list = c(
            "mr_ivw", "mr_egger_regression", "mr_weighted_median"
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

## heterogeneity test
mr_hete <- lapply(
  mr_data,
  function(x) {
    lapply(
      x,
      function(y) {
        mr_heterogeneity(y) %>%
          transform(
            Q_f = sprintf("%.2f", Q),
            p_f = ifelse(
              Q_pval < 0.001,
              sprintf("%.3e", Q_pval),
              sprintf("%.3f", Q_pval)
            )
          )
      }
    )
  }
)

## pleiotropy test by MR-Egger intercept
mr_plei <- lapply(
  mr_data,
  function(x) {
    lapply(
      x,
      function(y) {
        mr_pleiotropy_test(y) %>%
          transform(
            egger_intercept_f = sprintf("%.5f", egger_intercept),
            se_f = sprintf("%.5f", se),
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

## pleiotropy test by MR-PRESSO
cl <- makeCluster(detectCores())
mr_presso_raw <- lapply(
  mr_data,
  function(x) {
    parLapply(
      cl, x,
      function(y) {
        library(TwoSampleMR)
        set.seed(1)
        run_mr_presso(y)
      }
    )
  }
)
stopCluster(cl)

mr_presso <- lapply(
  mr_presso_raw,
  function(x) {
    lapply(
      x,
      function(y) {
        res_g <- y[[1]][["MR-PRESSO results"]][["Global Test"]]
        res_d <- y[[1]][["MR-PRESSO results"]][["Distortion Test"]]
        Global_Test_RSSobs <- sprintf("%.2f", res_g[["RSSobs"]])
        Global_Test_p <- sprintf("%.3f", res_g[["Pvalue"]])
        Outlier_No <- ifelse(
          is.null(res_d),
          NA,
          switch(class(res_d[["Outliers Indices"]]),
            "integer" = length(res_d[["Outliers Indices"]]),
            paste(res_d[["Outliers Indices"]], collapse = "|")
          )
        )
        Distortion_Coefficient <- ifelse(
          is.null(res_d),
          NA,
          sprintf("%.2f", res_d[["Distortion Coefficient"]])
        )
        Distortion_Test_p <- ifelse(
          is.null(res_d),
          NA,
          sprintf("%.3f", res_d[["Pvalue"]])
        )
        data.frame(
          Global_Test_RSSobs,
          Global_Test_p,
          Outlier_No,
          Distortion_Coefficient,
          Distortion_Test_p
        )
      }
    )
  }
)

## leave-one-out plots
mr_loo_plot <- lapply(
  mr_data,
  function(x) {
    lapply(
      x,
      function(y) {
        mr_leaveoneout(y) %>% geom_loo()
      }
    )
  }
)

## combine leave-one-out plots of pan-cancer
loo_pan <- mr_loo_plot[["OS"]][["All_sites"]] +
  mr_loo_plot[["CSS"]][["All_sites"]] +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

# data saving ----

## save MR results
write.xlsx(
  lapply(
    mr_res,
    function(x) {
      df <- data.frame()
      for (i in names(x)) {
        df <- rbind(df, x[[i]])
      }
      df
    }
  ),
  "10/MR_results.xlsx",
  TRUE
)

## save heterogeneity test results
write.xlsx(
  lapply(
    mr_hete,
    function(x) {
      df <- data.frame()
      for (i in names(x)) {
        df <- rbind(df, x[[i]])
      }
      df
    }
  ),
  "10/MR_heterogeneity.xlsx",
  TRUE
)

## save MR-Egger intercept test results
write.xlsx(
  lapply(
    mr_plei,
    function(x) {
      df <- data.frame()
      for (i in names(x)) {
        df <- rbind(df, x[[i]])
      }
      df
    }
  ),
  "10/MR_pleiotropy.xlsx",
  TRUE
)

## save MR-PRESSO results
write.xlsx(
  lapply(
    mr_presso,
    function(x) {
      df <- data.frame()
      for (i in names(x)) {
        dfi <- transform(x[[i]], cancer_type = cancer_names[[i]])
        df <- rbind(df, dfi)
      }
      df
    }
  ),
  "10/MR_presso.xlsx",
  TRUE
)

# plots saving ----

ggsave(
  "10/loo_pan-cacner.pdf",
  loo_pan,
  height = 7,
  width = 14,
  device = "pdf"
)
