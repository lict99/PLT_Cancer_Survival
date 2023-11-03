# env settings ----

rm(list = ls())
gc()

options(verbose = FALSE)
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

library(magrittr)
library(progress)

# run R scripts ----

r_files <- list.files(pattern = "^\\d\\d.+\\.R$") %>%
  sort()

pb <- progress_bar$new(
  format = "Running [:bar] :percent in :elapsed",
  total = length(r_files),
  clear = FALSE
)

one_step <- new.env(parent = globalenv())

cat("Start running ...\n")
for (fs in seq_along(r_files)) {
  pb$tick()
  source(r_files[fs], local = one_step, echo = FALSE, encoding = "UTF-8") %>%
    suppressMessages() %>%
    suppressWarnings() %>%
    capture.output() %>%
    invisible()
  rm(list = ls(envir = one_step), envir = one_step)
  gc()
}
cat("Mission complete!\n")

# TODO: add rm() and gc()
# TODO: add src/caches directory
