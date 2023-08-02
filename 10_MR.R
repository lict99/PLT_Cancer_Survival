# env settings ----

library(magrittr)
library(TwoSampleMR)
library(parallel)
library(openxlsx)

load("00/cancer_names.RData")
load("09/MR_harmonised_data.RData")

dir.create("10", FALSE)

set.seed(1)
cl <- makeCluster(detectCores())

# MR ----

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

mr_loo_plot <- lapply(
  mr_data,
  function(x) {
    lapply(
      x,
      function(y) {
        z <- mr_leaveoneout(y) %>% mr_leaveoneout_plot()
        z
      }
    )
  }
)

# data saving ----

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
