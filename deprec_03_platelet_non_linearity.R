
# env settings ------------------------------------------------------------

library(magrittr)
library(rms)
library(ggplot2)
library(cowplot)
load("00/ICD_of_cancers.RData")
load("00/ICD_of_cancers_excl_bood.RData")
load("01/tidy_data_diagnosis_after_attending.RData")
dir.create("03", FALSE)

# non linearity plot (rcs) ----------------------------------------

## all cancer type
plot_data <- list()
for (i in c(0, 365.25 / 2, 365.25)) {
  data_select <- tidy_data_dia_after_att[
    (tidy_data_dia_after_att$Date_of_diagnosis - tidy_data_dia_after_att$Date_attending) >= i,
    c("Date_of_diagnosis", "Date_end", "OS", "Age_at_recruitment", "Sex", "asprin", "SMOKING_STATUS", "ALCOHOL_STATUS", "BMI", "Townsend_deprivation_index.TDI._at_recruitment", "Race", "ICD", "Date_attending", "Platelet")
  ]
  y <- Surv(
    time = data_select$Date_end - data_select$Date_of_diagnosis,
    event = data_select$OS == 1
  )
  for (n in seq_along(ICD_rexp)) {
    cph.fit <- cph(
      y ~ rcs(Platelet, 3) + Age_at_recruitment + Sex + asprin + SMOKING_STATUS + ALCOHOL_STATUS + BMI + Townsend_deprivation_index.TDI._at_recruitment + Race,
      data_select,
      subset = grepl(ICD_rexp[[n]], data_select$ICD),
      x = T,
      y = T,
      surv = T
    )
    d <- datadist(data_select)
    options(datadist = "d")
    try(
      {
        effect <- Predict(cph.fit, Platelet, fun = exp, ref.zero = T)
        plot_data[[paste0(i, "_", n)]] <- effect
      },
      silent = T
    )
  }
}

plot_names <- names(plot_data)
plot_geom <- list()
for (j in plot_names) {
  p <- ggplot() +
    geom_ribbon(aes(x = Platelet, ymin = lower, ymax = upper), plot_data[[j]], fill = "grey", alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_line(aes(Platelet, yhat), plot_data[[j]], linewidth = 1) +
    xlab(expression(paste("Platelet Cutoff (", 10^9, "/L)", sep = ""))) +
    ylab("HR (95%CI)") +
    labs(subtitle = paste0(
      "Cancer Type: ", names(ICD_rexp)[as.numeric(strsplit(j, "_")[[1]][2])],
      "\n",
      "Lag Time: ", round(as.numeric(strsplit(j, "_")[[1]][1])), " day(s)"
    )) +
    theme_classic()
  plot_geom[[j]] <- p
}
p_merge <- plot_grid(plotlist = plot_geom, align = "hv")
save_plot("03/non_linear_rcs.pdf", p_merge, base_height = 18, base_width = 28)

## excl. blood cancer
plot_data <- list()
for (i in c(0, 365.25/2, 365.25)) {
  data_select <- tidy_data_dia_after_att[
    (tidy_data_dia_after_att$Date_of_diagnosis - tidy_data_dia_after_att$Date_attending) >= i,
    c("Date_of_diagnosis", "Date_end", "OS", "Age_at_recruitment", "Sex", "asprin", "SMOKING_STATUS", "ALCOHOL_STATUS", "BMI", "Townsend_deprivation_index.TDI._at_recruitment", "Race", "ICD", "Date_attending", "Platelet")
  ]
  y <- Surv(
    time = data_select$Date_end - data_select$Date_of_diagnosis,
    event = data_select$OS == 1
  )
  for (n in seq_along(ICD_rexp_excl_bood)) {
    cph.fit <- cph(
      y ~ rcs(Platelet, 3) + Age_at_recruitment + Sex + asprin + SMOKING_STATUS + ALCOHOL_STATUS + BMI + Townsend_deprivation_index.TDI._at_recruitment + Race,
      data_select,
      subset = !(grepl(ICD_rexp_excl_bood[[n]], data_select$ICD)),
      x = T,
      y = T,
      surv = T
    )
    d <- datadist(data_select)
    options(datadist = "d")
    try(
      {
        effect <- Predict(cph.fit, Platelet, fun = exp, ref.zero = T)
        plot_data[[paste0(i, "_", n)]] <- effect
      },
      silent = T
    )
  }
}

plot_names <- names(plot_data)
plot_geom <- list()
for (j in plot_names) {
  p <- ggplot() +
    geom_ribbon(aes(x = Platelet, ymin = lower, ymax = upper), plot_data[[j]], fill = "grey", alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_line(aes(Platelet, yhat), plot_data[[j]], linewidth = 1) +
    xlab(expression(paste("Platelet Cutoff (", 10^9, "/L)", sep = ""))) +
    ylab("HR (95%CI)") +
    labs(subtitle = paste0(
      "Cancer Type: ", names(ICD_rexp_excl_bood)[as.numeric(strsplit(j, "_")[[1]][2])],
      "\n",
      "Lag Time: ", round(as.numeric(strsplit(j, "_")[[1]][1])), " day(s)"
    )) +
    theme_classic()
  plot_geom[[j]] <- p
}
p_merge <- plot_grid(plotlist = plot_geom, align = "hv")
save_plot("03/non_linear_rcs_excl_blood.pdf", p_merge, base_height = 6, base_width = 8)
