# survival outcomes defined by ICD codes ----

extract_survival <- function(
    death_codes,
    diag_codes = cancer_ICD_codes,
    death_df = UKb_death,
    diag_df = UKb_diagnosis,
    when_end = date_end,
    baseline = UKb_baseline) {
  if (!identical(names(diag_codes), names(death_codes))) {
    stop("Invalid diag_codes and death_codes!")
  }
  # extract diagnosis
  diag_list <- lapply(
    diag_codes,
    function(x) {
      subset(diag_df, grepl(x, ICD_diagnosis))
    }
  )
  # extract events
  death_list <- lapply(
    death_codes,
    function(x) {
      subset(death_df, grepl(x, ICD10_death))
    }
  )
  # extract death date
  when_end_by_death <- tapply(
    death_df$date_death %>% as.Date(),
    death_df$eid,
    min
  ) %>%
    data.frame(last_fu = ., eid = names(.)) %>%
    transform(last_fu = as.Date(last_fu))
  # merge diagnosis and events
  list <- list()
  for (i in names(diag_codes)) {
    diag <- aggregate(
      diag_list[[i]][, c("ICD_diagnosis", "date_diagnosis")],
      list(eid = diag_list[[i]][, "eid"]),
      function(x) {
        if (class(x) == "character") {
          unique(x) %>% paste(collapse = " | ")
        } else if (class(x) == "Date") {
          min(x)
        } else {
          stop("Invalid data class!")
        }
      }
    ) %>%
      transform(date_diagnosis = as.Date(date_diagnosis))
    death <- aggregate(
      death_list[[i]][, c("ICD10_death", "date_death")],
      list(eid = death_list[[i]][, "eid"]),
      function(x) {
        if (class(x) == "character") {
          unique(x) %>% paste(collapse = " | ")
        } else if (class(x) == "Date") {
          min(x)
        } else {
          stop("Invalid data class!")
        }
      }
    ) %>%
      transform(date_death = as.Date(date_death))
    list[[i]] <- merge(diag, death, by = "eid", all.x = TRUE) %>%
      merge(when_end_by_death, by = "eid", all.x = TRUE) %>%
      transform(
        fu_event = ifelse(is.na(date_death), 0, 1),
        last_fu = ifelse(
          is.na(last_fu),
          as.character(when_end),
          as.character(last_fu)
        ) %>% as.Date()
      ) %>%
      transform(
        fu_time = difftime(
          last_fu,
          date_diagnosis,
          units = "days"
        )
      ) %>%
      merge(baseline, by = "eid", all.x = TRUE) %>%
      transform(
        lag_time = difftime(
          date_diagnosis,
          date_attending,
          units = "days"
        ),
        platelet_per_100 = platelet / 100,
        platelet_300 = ifelse(platelet >= 300, "YES", "NO") %>% factor(),
        platelet_400 = ifelse(platelet >= 400, "YES", "NO") %>% factor()
      ) %>%
      subset(fu_time >= 0) # delete cases with invalid follow-up time
  }
  list
}
