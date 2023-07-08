# env settings ----

library(magrittr)
library(table1)
library(openxlsx)

load(file = "01/whole_cancer_data_for_Cox.RData")
load(file = "00/cancer_ICD_codes_with_attr.RData")
source("functions/Cox_regression.R")

dir.create("02", FALSE)

# statistics ----

vars_num <- c("platelet", "age", "bmi", "TDI")
vars_fct <- c(
  "platelet_300",
  "platelet_400",
  "sex",
  "aspirin",
  "smoking_status",
  "alcohol_status",
  "race"
)

## compute manually
bsl <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = c(
        vars_num,
        vars_fct,
        "fu_event",
        "lag_time"
      ),
      lagtime = c(-Inf, Inf)
    ) %>%
      lapply(
        function(y) {
          res <- list()
          res[["stats"]] <- by(
            y, ~ lag_time >= 0,
            function(z) {
              list <- list()
              for (i in colnames(z)) {
                if (i %in% vars_num) {
                  list[[i]] <- list(
                    desc = paste0(
                      median(z[, i], na.rm = TRUE) %>% sprintf("%.3f", .),
                      "(",
                      quantile(z[, i], 0.25, na.rm = TRUE) %>%
                        sprintf("%.3f", .),
                      ", ",
                      quantile(z[, i], 0.75, na.rm = TRUE) %>%
                        sprintf("%.3f", .),
                      ")"
                    ),
                    na = paste0(
                      sum(is.na(z[, i])),
                      "(",
                      (sum(is.na(z[, i])) / length(z[, i])) %>%
                        sprintf("%.1f", .),
                      ")"
                    )
                  )
                } else if (i %in% vars_fct) {
                  list[[i]] <- list(
                    table = table(z[, i], useNA = "ifany"),
                    prop = prop.table(table(z[, i], useNA = "ifany")) %>%
                      multiply_by(100) %>%
                      sprintf("%.1f", .)
                  )
                } else if (i %in% "fu_event") {
                  list[[i]] <- list(
                    nevent = sum(z[, i]),
                    n = length(z[, i])
                  )
                }
              }
              list
            },
            simplify = FALSE
          )
          res[["test"]] <- (function() {
            list <- list()
            for (j in colnames(y)) {
              if (j %in% vars_num) {
                list[[j]] <- wilcox.test(
                  as.formula(paste(j, "~ lag_time >= 0")),
                  data = y,
                  paired = FALSE
                )[["p.value"]]
              } else if (j %in% vars_fct) {
                list[[j]] <- chisq.test(
                  table(
                    transform(
                      extract(y, c(j, "lag_time")),
                      lag_time = lag_time >= 0
                    ),
                    useNA = "no"
                  ),
                  correct = FALSE
                )[["p.value"]]
              }
            }
            list
          })()
          res
        }
      )
  }
)

## produce tables
tables <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = c(
        vars_num,
        vars_fct,
        "fu_event",
        "lag_time"
      ),
      lagtime = c(-Inf, Inf)
    ) %>%
      lapply(
        function(y) {
          transform(
            y,
            lag_time = ifelse(
              lag_time >= 0,
              "lag_time >= 0",
              "lag_time < 0"
            ) %>%
              factor(),
            fu_event = factor(
              fu_event,
              levels = c(0, 1),
              labels = c("Censored events", "Failure events")
            )
          ) %>%
            table1(
              ~ . | lag_time,
              data = .,
              overall = FALSE
            ) %>%
            as.data.frame()
        }
      )
  }
)

# data saving ----

save(bsl, file = "02/baseline_stats_test.RData")
write.xlsx(tables[["OS"]], file = "02/table1_for_OS.xlsx")
write.xlsx(tables[["CSS"]], file = "02/table1_for_CSS.xlsx")
