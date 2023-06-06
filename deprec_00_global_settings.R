
# env settings ------------------------------------------------------------
library(magrittr)
dir.create("00", FALSE)

# common variables --------------------------------------------------------

## ICD codes
ICD_rexp <- list(
  pan_cancer = "^C[0-8][0-9]|^C9[0-7]|^1[4-9][0-9]|^20[0-8]",
  lung_cancer = "^C3[34]|^162",
  colorectal_cancer = "^C1[89]|^C20|^153[0-9]|^154[01]",
  liver_cancer = "^C22|^155",
  stomach_cancer = "^C16|^151",
  breast_cancer = "^C50|^17[45]",
  esophageal_cancer = "^C15|^150",
  pancreatic_cancer = "^C25|^157",
  prostate_cancer = "^C61|^185",
  cervical_cancer = "^C53|^180",
  leukemia = "^C9[1-5]|^2089",
  endometrial_cancer = "^C541|^1820",
  glioblastoma = "^C71|^191",
  lymphoma = "^C8[1-9]|^C9[0-6]|^2028",
  ovarian_cancer = "^C56|^183|^220",
  bladder_cancer = "^C67|^188"
) %>% 
  lapply(X = ., FUN = `attr<-`, which = "incl", value = TRUE)

## ICD codes used for excluding some cancers
ICD_rexp_excl_cancers <- list(
  excl_blood = "^C8[1-9]|^C9[0-6]|^2028|^2089"
) %>% 
  lapply(X = ., FUN = `attr<-`, which = "incl", value = FALSE)

## ICD codes combination
ICD_with_attr <- c(ICD_rexp_excl_cancers, ICD_rexp)

# functions ---------------------------------------------------------------

## Cox regression loop
UKb_Cox_loop <- function(data,
                         start = 0,
                         end,
                         event,
                         eventcode = 1,
                         covars,
                         ICD,
                         ICDrexp,
                         attending,
                         lag = 0,
                         target) {
  requireNamespace("survival")
  
  cox.df <- data.frame()
  
  for (i in lag) {
    data_lag <- subset(data, (data[,start] - data[,attending]) >= i)
    
    y <- survival::Surv(
      time = data_lag[,end] - data_lag[,start],
      event = data_lag[,event] == eventcode
    )
    
    for (j in target) {
      if (!is.null(covars)) {
        formula <- as.formula(
          paste0(
            "y ~ ",
            j,
            " + ",
            paste(covars, collapse = " + ")
          )
        )
      } else {
        formula <- as.formula(paste0("y ~ ", j))
      }
      
      for (n in seq_along(ICDrexp)) {
        subset <- grepl(ICDrexp[[n]], data_lag[,ICD])
        
        cox.fit <- survival::coxph(
          formula = formula,
          data = data_lag,
          subset = subset
        )
        
        summary.cox <- summary(cox.fit)
        
        target.ci <- grep(j, rownames(summary.cox[["conf.int"]]))
        target.coef <- grep(j, rownames(summary.cox[["coefficients"]]))
        
        use_data <- na.omit(
          data_lag[subset,c(start, end, event, covars, ICD, attending, j)],
        )
        
        n_target <- if (
          identical(
            nrow(use_data), 
            cox.fit[["n"]]
          )
        ) {
          if (class(data_lag[,j]) == "factor") {
            length(which(use_data[,j] == "YES"))
          } else {NA}
        } else {
          stop("check your data!")
        }
        
        target_df <- data.frame(
          cancer_type = names(ICDrexp)[n],
          ICD = paste(
            sort(unique(grep(ICDrexp[[n]], data_lag[,ICD], value = T))), 
            collapse = "|"
          ),
          n = cox.fit[["n"]],
          n_target = n_target,
          percent_target = n_target/(cox.fit[["n"]]),
          nevent = cox.fit[["nevent"]],
          missing = length(cox.fit[["na.action"]]),
          lag_time = i,
          covars = paste(covars, collapse = "|"),
          target = rownames(summary.cox[["conf.int"]])[target.ci],
          HR = summary.cox[["conf.int"]][target.ci,"exp(coef)"],
          lower.95 = summary.cox[["conf.int"]][target.ci,"lower .95"],
          upper.95 = summary.cox[["conf.int"]][target.ci,"upper .95"],
          coef = summary.cox[["coefficients"]][target.coef,"coef"],
          P = summary.cox[["coefficients"]][target.coef,"Pr(>|z|)"]
        )
        
        cox.df <- rbind(cox.df, target_df, make.row.names = F)
      }
    }
  }
  cox.df <- transform(cox.df, fdr = p.adjust(P, method = "fdr"))
  return(cox.df)
}
save(UKb_Cox_loop, file = "00/fun_UKb_Cox_loop.RData")

## Cox regression loop excluded some cancer type
UKb_Cox_loop_excluded <- function(
    data,
    start = 0,
    end,
    event,
    eventcode = 1,
    covars,
    ICD,
    ICDrexp,
    attending,
    lag = 0,
    target) {
  requireNamespace("survival")
  
  cox.df <- data.frame()
  
  for (i in lag) {
    data_lag <- subset(data, (data[,start] - data[,attending]) >= i)
    
    y <- survival::Surv(
      time = data_lag[,end] - data_lag[,start],
      event = data_lag[,event] == eventcode
    )
    
    for (j in target) {
      if (!is.null(covars)) {
        formula <- as.formula(
          paste0(
            "y ~ ",
            j,
            " + ",
            paste(covars, collapse = " + ")
          )
        )
      } else {
        formula <- as.formula(paste0("y ~ ", j))
      }
      
      for (n in seq_along(ICDrexp)) {
        subset <- !(grepl(ICDrexp[[n]], data_lag[,ICD]))
        
        cox.fit <- survival::coxph(
          formula = formula,
          data = data_lag,
          subset = subset
        )
        
        summary.cox <- summary(cox.fit)
        
        target.ci <- grep(j, rownames(summary.cox[["conf.int"]]))
        target.coef <- grep(j, rownames(summary.cox[["coefficients"]]))
        
        use_data <- na.omit(
          data_lag[subset,c(start, end, event, covars, ICD, attending, j)],
        )
        
        n_target <- if (
          identical(
            nrow(use_data), 
            cox.fit[["n"]]
          )
        ) {
          if (class(data_lag[,j]) == "factor") {
            length(which(use_data[,j] == "YES"))
          } else {NA}
        } else {
          stop("check your data!")
        }
        
        target_df <- data.frame(
          cancer_type = names(ICDrexp)[n],
          ICD = paste(
            sort(unique(grep(ICDrexp[[n]], data_lag[,ICD], value = T, invert = T))), 
            collapse = "|"
          ),
          n = cox.fit[["n"]],
          n_target = n_target,
          percent_target = n_target/(cox.fit[["n"]]),
          nevent = cox.fit[["nevent"]],
          missing = length(cox.fit[["na.action"]]),
          lag_time = i,
          covars = paste(covars, collapse = "|"),
          target = rownames(summary.cox[["conf.int"]])[target.ci],
          HR = summary.cox[["conf.int"]][target.ci,"exp(coef)"],
          lower.95 = summary.cox[["conf.int"]][target.ci,"lower .95"],
          upper.95 = summary.cox[["conf.int"]][target.ci,"upper .95"],
          coef = summary.cox[["coefficients"]][target.coef,"coef"],
          P = summary.cox[["coefficients"]][target.coef,"Pr(>|z|)"]
        )
        
        cox.df <- rbind(cox.df, target_df, make.row.names = F)
      }
    }
  }
  cox.df <- transform(cox.df, fdr = p.adjust(P, method = "fdr"))
  return(cox.df)
}
save(UKb_Cox_loop_excluded, file = "00/fun_UKb_Cox_loop_excluded.RData")

# data saving -------------------------------------------------------------

save(ICD_with_attr, file = "00/ICD_with_attr.RData")
