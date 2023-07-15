# env settings ----

library(magrittr)
library(table1)
library(openxlsx)

load(file = "01/whole_cancer_data_for_Cox.RData")
load(file = "00/cancer_ICD_codes_with_attr.RData")
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
              overall = FALSE,
              extra.col = NULL
            ) %>%
            as.data.frame()
        }
      )
  }
)

# data saving ----

write.xlsx(tables[["OS"]], file = "02/table1_for_OS.xlsx")
write.xlsx(tables[["CSS"]], file = "02/table1_for_CSS.xlsx")
