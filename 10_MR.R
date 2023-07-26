# env settings ----

library(magrittr)
library(TwoSampleMR)

load("09/MR_harmonised_data.RData")

dir.create("10", FALSE)

set.seed(1)

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
    lapply(
      x,
      function(y) {
        run_mr_presso(y)
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

save(mr_res, file = "10/MR_results.RData")
write.xlsx(mr_res[["OS"]], "10/MR_results_OS.xlsx", TRUE)
write.xlsx(mr_res[["CSS"]], "10/MR_results_CSS.xlsx", TRUE)
write.xlsx(mr_hete[["OS"]], "10/MR_heterogeneity_OS.xlsx", TRUE)
write.xlsx(mr_hete[["CSS"]], "10/MR_heterogeneity_CSS.xlsx", TRUE)
write.xlsx(mr_plei[["OS"]], "10/MR_pleiotropy_OS.xlsx", TRUE)
write.xlsx(mr_plei[["CSS"]], "10/MR_pleiotropy_CSS.xlsx", TRUE)
