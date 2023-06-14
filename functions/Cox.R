# env settings ----

library(magrittr)
library(survival)

# functions used in Cox proportional hazard regression ----

## extract data by vars and filter by "sex"
extract_Cox_data <- function(
    data_list = whole_cancer_data,
    vars,
    lagtime_col = "lag_time",
    lagtime = NULL,
    ICD_codes = cancer_ICD_codes) {
  if (!is.numeric(lagtime) && length(lagtime) != 2) {
    stop("Lag time window is not two numbers!")
  }
  list <- lapply(
    data_list,
    function(x) {
      x[x[, lagtime_col] >= min(lagtime) & x[, lagtime_col] < max(lagtime), ]
    }
  ) %>%
    lapply(
      function(x) {
        x[, unique(c(vars, "sex"))]
      }
    )
  for (i in names(ICD_codes)) {
    if (length(attr(ICD_codes[[i]], "sex") != 0)) {
      if (attr(ICD_codes[[i]], "sex") == "male") {
        list[[i]] <- (list[[i]][, "sex"] == "Male") %>%
          extract(list[[i]], ., vars[vars != "sex"], drop = FALSE)
      } else if (attr(ICD_codes[[i]], "sex") == "Female") {
        list[[i]] <- (list[[i]][, "sex"] == "Female") %>%
          extract(list[[i]], ., vars[vars != "sex"], drop = FALSE)
      }
    }
  }
  list
}

## run Cox regression per lag time, data is from extract_Cox_data() call
run_Cox_per_lag <- function(
    ext_data,
    futime_col = "fu_time",
    event_col = "cancer_death",
    eventcode = 1) {
  fit <- lapply(
    ext_data,
    function(x) {
      f <- as.formula(
        paste0(
          "Surv(", futime_col, ", ", event_col, "==", eventcode, ") ~ ."
        )
      )
      coxph(f, data = x, singular.ok = TRUE)
    }
  ) %>%
    lapply(
      function(x) {
        list(
          smr = summary(x),
          zph = try(cox.zph(x)[["table"]], silent = TRUE)
        )
      }
    )
  fit
}

## extract summary table from run_Cox_per_lag() call
extract_smr_data <- function(
    run_data,
    target = "platelet") {
  list <- lapply(
    run_data,
    function(x) {
      zph_p <- x[["zph"]]["GLOBAL", "p"]
      smr <- x[["smr"]]
      df <- data.frame()
    }
  )
}

## run Cox regression for multiple lag time values
run_Cox_loop <- function(mlagtime) {
  extract_Cox_data() %>%
    run_Cox_per_lag() %>%
    extract_smr_data()
}
