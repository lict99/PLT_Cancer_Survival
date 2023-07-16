# env settings ----

library(magrittr)

# run R scripts ----

r_files <- list.files(pattern = "^0[0-6].+R$") %>%
  sort()

for (fs in r_files) {
  message("Start running \"", fs, "\".")
  source(fs, echo = FALSE, encoding = "UTF-8")
  message("Successfully run \"", fs, "\".")
  rm(list = ls())
  gc()
}

