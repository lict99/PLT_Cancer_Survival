# functions used in Cox proportional harzard regression ----

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

## run Cox by multiple lag time in each cancer
run_Cox_loop <- function(
    data_list = whole_cancer_data, # nested extract_Cox_data() call
    futime_col = "fu_time",
    event_col = "cancer_death",
    eventcode = 1,
    vars,
    lagtime_col = "lag_time",
    mlagtime,
    ICD_codes = cancer_ICD_codes) {
  run_Cox_per_cancer <- function(lagtime, ...) {
    list <- extract_Cox_data(
      data_list = data_list,
      vars = c(futime_col, event_col, vars),
      lagtime_col = lagtime_col,
      lagtime = lagtime,
      ICD_codes = ICD_codes
    )
    fit <- lapply(
      list,
      function(x) {
        f <- as.formula(
          paste0(
            "Surv(", futime_col, ", ", event_col, "==", eventcode, ") ~ ."
          )
        )
        coxph(f, data = x, singular.ok = TRUE)
      }
    )
    lapply(
      fit,
      function(x) {
        list(
          smr = summary(x),
          zph = try(cox.zph(x)[["table"]], silent = TRUE)
        )
      }
    )
  }
  res <- list()
  for (i in seq_along(mlagtime)) {
    res[[paste("lagtime", min(mlagtime[[i]]), "days", sep = "_")]] <- do.call(
      run_Cox_per_cancer,
      list(
        lagtime = mlagtime[[i]],
        data_list = data_list,
        futime_col = futime_col,
        event_col = event_col,
        eventcode = eventcode,
        vars = vars,
        lagtime_col = lagtime_col,
        ICD_codes = ICD_codes
      )
    )
  }
  res
}

extract_loop_data <- function() {

}
