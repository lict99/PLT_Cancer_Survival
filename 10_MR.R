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
        mr_heterogeneity(y)
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
        mr_pleiotropy_test(y)
      }
    )
  }
)

mr_presso <- lapply(
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
