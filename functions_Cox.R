# functions used in Cox proportional harzard regression ----

## extract data by vars and filter by "sex"
extract_Cox_data <- function(
    data_list = whole_cancer_data,
    vars,
    lagtime_col = "lag_time",
    lagtime = NULL,
    ICD_codes = cancer_ICD_codes) {
  if (!is.numeric(lagtime) && length(lagtime) != 2) {
    stop("Lag time window is not two numbers")
  }
  list <- lapply(
    data_list,
    function(x) {
      x[x[, lagtime_col] >= min(lagtime) & x[, lagtime_col] < max(lagtime), ]
    }
  ) %>%
    lapply(
      function(x) {
        x[, vars]
      }
    )
  for (i in names(ICD_codes)) {
    if (length(attr(ICD_codes[[i]], "sex") != 0)) {
      if (attr(ICD_codes[[i]], "sex") == "male") {
        list[[i]] <- (list[[i]][, "sex"] == "Male") %>%
          extract(list[[i]], ., vars[vars != "sex"])
      } else if (attr(ICD_codes[[i]], "sex") == "Female") {
        list[[i]] <- (list[[i]][, "sex"] == "Female") %>%
          extract(list[[i]], ., vars[vars != "sex"])
      }
    }
  }
  list
}

## run Cox by multiple lag time in each cancer
run_Cox_loop <- function() {
  run_Cox_per_cancer <- function() {

  }
}
