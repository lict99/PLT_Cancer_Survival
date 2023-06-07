# env settings ------------------------------------------------------------

library(magrittr)

dir.create("00", FALSE)

# common variables --------------------------------------------------------

## end of follow-up time
date_end <- as.Date("2021-07-01")

## cancer ICD codes
ICD_incl <- list(
  All_sites = "^C[0-8][0-9]|^C9[0-7]|^1[4-9][0-9]|^20[0-8]",
  Female_breast = "^C50|^174" %>%
    set_attr("sex", "female"),
  Lung = "^C3[34]|^162",
  Colorectum = "^C1[89]|^C2[01]|^15[34]",
  Prostate = "^C61|^185" %>%
    set_attr("sex", "male"),
  Nonmelanoma_of_skin = "^C44|^173",
  Stomach = "^C16|^151",
  Liver = "^C22|^155",
  Cervix_uteri = "^C53|^180" %>%
    set_attr("sex", "female"),
  Esophagus = "^C15|^150",
  Thyroid = "^C73|^193",
  Bladder = "^C67|^188",
  Non_Hodgkin_lymphoma = "^C8[2-6]|^C96|^20[02]",
  Pancreas = "^C25|^157",
  Leukemia = "^C9[1-5]|^20[4-8]",
  Kidney = "^C6[45]|^189",
  Corpus_uteri = "^C54|^182" %>%
    set_attr("sex", "female"),
  Lip_oral_cavity = "^C0[0-6]|^14[01345]",
  Melanoma_of_skin = "^C43|^172",
  Ovary = "^C56|^1830" %>%
    set_attr("sex", "female"),
  Brain_nervous_system = "^C7[0-2]|^19[12]",
  Larynx = "^C32|^161",
  Multiple_myeloma = "^C88|^C90|^203",
  Nasopharynx = "^C11|^147",
  Gallbladder = "^C23|^1560",
  Oropharynx = "^C09|^C10|^146",
  Hypopharynx = "^C1[23]|^148",
  Hodgkin_lymphoma = "^C81|^201",
  Testis = "^C62|^186" %>%
    set_attr("sex", "male"),
  Salivary_glands = "^C0[78]|^142",
  # Anus = "", # merged to colorectum
  Vulva = "^C51|^1844" %>%
    set_attr("sex", "female"),
  Penis = "^C60|^187" %>%
    set_attr("sex", "male"),
  Kaposi_sarcoma = "^C46|^176",
  Mesothelioma = "^C45", # no ICD-9 code exists
  Vagina = "^C52|^1840" %>%
    set_attr("sex", "female")
) %>%
  lapply(X = ., FUN = `attr<-`, which = "incl", value = TRUE)

## cancer ICD codes excluding some types of cancers
ICD_excl <- list(
  E_nonmelanoma_skin = "^C44|^173",
  E_lymphoid_haematopoietic = "^C8[1-9]|^C9[0-6]|^20[0-8]"
) %>%
  lapply(X = ., FUN = `attr<-`, which = "incl", value = FALSE)

## combination of all ICD codes
cancer_ICD_codes <- c(ICD_incl, ICD_excl)

# functions ---------------------------------------------------------------

extract_Cox_data <- function(data_list = whole_cancer_data,
                             vars,
                             lagtime_col = "lag_time",
                             lagtime = -Inf,
                             ICD_codes = cancer_ICD_codes) {
  if (is.numeric(lagtime)) {
    data_list <- lapply(
      data_list,
      function(x) {
        x[x[, lagtime_col] >= lagtime, ]
      }
    )
  } else {
    stop("Lag time is not numeric!")
  }

  list <- lapply(
    data_list,
    function(x) {
      x[, vars]
    }
  )

  if (any(grepl("sex", vars, ignore.case = TRUE), na.rm = TRUE)) {
    for (i in names(cancer_ICD_codes)) {
      if (length(attr(cancer_ICD_codes[[i]], "sex")) != 0) {
        if (attr(cancer_ICD_codes[[i]], "sex") == "male") {
          list[[i]] <- list[[i]] %>%
            subset(subset = sex == "Male") %>%
            subset(select = -sex)
        } else if (attr(cancer_ICD_codes[[i]], "sex") == "female") {
          list[[i]] <- list[[i]] %>%
            subset(subset = sex == "Female") %>%
            subset(select = -sex)
        }
      }
    }
  }

  return(list)
}

run_Cox_per_cancer <- function(data_list = whole_cancer_data,
                               futime_col = "fu_time",
                               event_col = "cancer_death",
                               eventcode = 1,
                               vars,
                               lagtime_col = "lag_time",
                               lagtime = -Inf,
                               ICD_codes = cancer_ICD_codes,
                               return = "coxph") {
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
      coxph(f, data = x, singular.ok = TRUE, model = TRUE, x = TRUE)
    }
  )

  if (return == "coxph") {
    return(fit)
  }

  zph <- list()
  for (i in names(fit)) {
    tryCatch(
      zph[[i]] <- cox.zph(fit[[i]]),
      error = function(e) {
        "Error!"
      }
    )
  }

  if (return == "cox.zph") {
    return(zph)
  } else {
    return("Invalid return!")
  }
}

run_Cox_loop <- function(mlagtime,
                         vars,
                         ...) {
  list <- list()
  for (i in mlagtime) {
    list[[paste("lagtime", i, "days", sep = "_")]] <- do.call(
      run_Cox_per_cancer,
      c(list(lagtime = i, vars = vars), list(...))
    )
  }
  return(list)
}

extract_loop_data <- function(loop_data, target) {
  ls <- list()
  for (i in names(loop_data)) {
    smr <- loop_data[[i]] %>%
      lapply(summary) %>%
      sapply(
        function(x) {
          data.frame(
            target = target,
            covars = ,
            icd_codes = ,
          )
        }
      )
  }

  ls[[i]] <- x
}

# data saving -------------------------------------------------------------

save(date_end, file = "00/end_of_follow-up_time.RData")
save(cancer_ICD_codes, file = "00/cancer_ICD_codes_with_attr.RData")
save(
  extract_Cox_data,
  run_Cox_per_cancer,
  run_Cox_loop,
  file = "00/functions.RData"
)
