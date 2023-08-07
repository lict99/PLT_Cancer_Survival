# env settings ----

library(magrittr)
library(table1)
library(openxlsx)

load(file = "01/whole_cancer_data_for_Cox.RData")
load(file = "00/cancer_ICD_codes_with_attr.RData")
load(file = "00/cancer_names.RData")
source("functions/Cox_regression.R")

dir.create("02", FALSE)

# statistics ----

vars <- c(
  "platelet", "age", "sex", "bmi", "TDI", "aspirin",
  "smoking_status", "alcohol_status", "race",
  "fu_time", "fu_event", "lag_time"
)

## produce tables
tables <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = vars,
      lagtime = c(-Inf, Inf)
    ) %>%
      lapply(
        function(y) {
          transform(
            y,
            aspirin = factor(aspirin, levels = c("YES", "NO")),
            race = factor(race, levels = c("British", "Others")),
            lag_time = ifelse(
              lag_time >= 0,
              "lag_time >= 0",
              "lag_time < 0"
            ) %>%
              factor(),
            fu_time = as.numeric(fu_time) %>% divide_by(365.25),
            fu_event = factor(
              fu_event,
              levels = c(1, 0),
              labels = c("Death", "Censoring")
            )
          ) %>%
            table1(
              ~ . | lag_time,
              data = .,
              render.continuous = function(x) {
                with(
                  stats.default(x),
                  c(
                    "",
                    `Median (IQR)` = sprintf("%.2f (%.2f-%.2f)", MEDIAN, Q1, Q3)
                  )
                )
              },
              overall = FALSE,
              extra.col = NULL
            ) %>%
            as.data.frame()
        }
      )
  }
)

## calculate number and proportion
## all participants with cancer including eligible and ineligible patients
nprop <- list(
  lag_no_limit = extract_Cox_data(
    data_list = whole_cancer_data[["OS"]],
    vars = c("eid"),
    lagtime = c(-Inf, Inf)
  ) %>%
    {
      n <- sapply(., nrow)
      n_max <- max(n)
      prop <- sapply(., function(x) nrow(x) / n_max)
      data.frame(
        cancer = unlist(cancer_names[names(.)]),
        N = n,
        proportion = paste0(sprintf("%.1f", prop * 100), "%")
      )
    },
  lag_0 = extract_Cox_data(
    data_list = whole_cancer_data[["OS"]],
    vars = c("eid"),
    lagtime = c(0, Inf)
  ) %>%
    {
      n <- sapply(., nrow)
      n_max <- max(n)
      prop <- sapply(., function(x) nrow(x) / n_max)
      data.frame(
        cancer = unlist(cancer_names[names(.)]),
        N = n,
        proportion = paste0(sprintf("%.1f", prop * 100), "%")
      )
    }
)

# data saving ----

write.xlsx(tables[["OS"]], file = "02/table1_for_OS.xlsx")
write.xlsx(tables[["CSS"]], file = "02/table1_for_CSS.xlsx")
write.xlsx(nprop, file = "02/number_and_proportion.xlsx")
