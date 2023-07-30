# env settings ----

library(magrittr)

set.seed(1)

# run R scripts ----

r_files <- list.files(pattern = "^\\d\\d.+\\.R$") %>%
  sort()

for (fs in r_files) {
  message("Start running \"", fs, "\".")
  source(fs, echo = FALSE, encoding = "UTF-8")
  message("Successfully run \"", fs, "\".")
  rm(list = ls())
  gc()
}
