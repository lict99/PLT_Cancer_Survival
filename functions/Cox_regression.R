# functions used in Cox proportional hazard regression ----

## extract data by vars and filter by "sex"
extract_Cox_data <- function(
    data_list, # a list nested with different diseases data
    vars,
    lagtime_col = "lag_time",
    lagtime = NULL,
    ICD_codes = cancer_ICD_codes) {
  if (!is.numeric(lagtime) || length(lagtime) != 2) {
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
      } else if (attr(ICD_codes[[i]], "sex") == "female") {
        list[[i]] <- (list[[i]][, "sex"] == "Female") %>%
          extract(list[[i]], ., vars[vars != "sex"], drop = FALSE)
      }
    }
  }
  if (!is.element("sex", vars)) {
    list <- lapply(
      list,
      function(x) {
        x[, colnames(x) != "sex", drop = FALSE]
      }
    )
  }
  list
}

## run Cox regression per lag time, data is from extract_Cox_data() call
run_Cox_per_lag <- function(
    ext_data,
    futime_col = "fu_time",
    event_col = "fu_event",
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
    target) {
  list <- lapply(
    run_data,
    function(x) {
      smr <- x[["smr"]]
      data.frame(
        n = smr[["n"]],
        nevent = smr[["nevent"]],
        covars = rownames(smr[["conf.int"]]) %>%
          extract(!grepl(target, .)) %>%
          paste(collapse = " | "),
        target = rownames(smr[["conf.int"]]) %>%
          extract(grepl(target, .)),
        HR = smr[["conf.int"]] %>%
          extract(grepl(target, rownames(.)), "exp(coef)"),
        lower.95 = smr[["conf.int"]] %>%
          extract(grepl(target, rownames(.)), "lower .95"),
        upper.95 = smr[["conf.int"]] %>%
          extract(grepl(target, rownames(.)), "upper .95"),
        p = smr[["coefficients"]] %>%
          extract(grepl(target, rownames(.)), "Pr(>|z|)"),
        zph_p_target = tryCatch(
          x[["zph"]] %>%
            extract(grepl(target, rownames(.)), "p"),
          error = function(...) NA
        ) %>%
          ifelse(length(.) == 1, ., NA)
      ) %>%
        transform(
          HR_f = paste0(
            sprintf("%.2f", HR),
            " (",
            sprintf("%.2f", lower.95), "-", sprintf("%.2f", upper.95),
            ")"
          ),
          p_f = ifelse(p < 0.001, sprintf("%.3e", p), sprintf("%.3f", p))
        )
    }
  )
  df <- data.frame()
  for (i in seq_along(list)) {
    dfi <- data.frame(cancer_site = names(list)[i]) %>%
      cbind(list[[i]])
    df <- rbind(df, dfi)
  }
  df
}

## run Cox regression for multiple lag time values
run_Cox_loop <- function(
    mlagtime,
    vars,
    target,
    data_list,
    lagtime_col = "lag_time",
    ICD_codes = cancer_ICD_codes,
    futime_col = "fu_time",
    event_col = "fu_event",
    eventcode = 1) {
  df <- data.frame()
  for (i in seq_along(mlagtime)) {
    lag <- paste(
      min(mlagtime[[i]]),
      "to",
      max(mlagtime[[i]]),
      "day(s)"
    )
    dfi <- extract_Cox_data(
      data_list,
      vars,
      lagtime_col,
      lagtime = mlagtime[[i]],
      ICD_codes
    ) %>%
      run_Cox_per_lag(
        futime_col,
        event_col,
        eventcode
      ) %>%
      extract_smr_data(target) %>%
      transform(lag_time = lag)
    df <- rbind(df, dfi)
  }
  df
}
