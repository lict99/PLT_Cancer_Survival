
# env settings ------------------------------------------------------------

library(magrittr)
library(smoothHR)
library(smoothHR)
library(Hmisc)
library(ggplot2)
library(ggsci)
library(scales)
library(cowplot)

load("01/whole_cancer_data_for_Cox.RData")
load("00/functions.RData")
load("00/cancer_ICD_codes_with_attr.RData")

dir.create("04", FALSE)


# function definition -----------------------------------------------------

cal_hr <- function(data) {
  if (any(grepl("sex", colnames(data)))) {
    fit <- coxph(
      formula = Surv(fu_time, cancer_death == 1) ~ 
        rcspline.eval(platelet, nk = 3, inclx = TRUE) + sex +
        age + asprin + smoking_status + alcohol_status + bmi + TDI + race,
      data = data,
      x = TRUE,
      model = TRUE,
      singular.ok = TRUE
    )
  } else {
    fit <- coxph(
      formula = Surv(fu_time, cancer_death == 1) ~ 
        rcspline.eval(platelet, nk = 3, inclx = TRUE) + 
        age + asprin + smoking_status + alcohol_status + bmi + TDI + race,
      data = data,
      x = TRUE,
      model = TRUE,
      singular.ok = TRUE
    )
  }
  
  sumr <- summary(fit)
  
  hr <- smoothHR(data = data, coxfit = fit)
  
  hr_point <- predict.HR(
    hr,
    predictor = "platelet",
    prob = 0.5,
    # pred.value = NULL,
    prediction.values = seq(
      min(hr$dataset$platelet), 
      max(hr$dataset$platelet), 
      length.out = 100
    ),
    conf.level = 0.95
  ) %>%
    as.data.frame() %>% 
    within(
      {
        HR <- exp(LnHR)
        lower_95 <- exp(`lower .95`)
        upper_95 <- exp(`upper .95`)
      }
    )
  
  list(dataset = hr$dataset, sumr = sumr, hr_point = hr_point)
}

geom_hr <- function(data) {
  colors <- pal_nejm(palette = c("default"), alpha = 1)(8) %T>% show_col() 
  dataset <- data$dataset
  p <- data$sumr$logtest[3]
  point <- data$hr_point
  dst <- density(dataset$platelet)
  dst_lmt <- c(0, max(dst$y))
  hr_lmt <- c(
    0, 
    ifelse(
      (max(point$upper_95)/max(point$HR)) > 2, 
      max(point$HR), 
      max(point$upper_95)
    )
  )
  dst_y_rsc <- rescale(dst$y, hr_lmt, dst_lmt)
  
  col_dst <- colors[5]
  col_hr <- colors[1]
  col_ci <- "black"
  
  p1 <- ggplot() +
    geom_ribbon(
      aes(x = dst$x, ymin = min(dst_y_rsc), ymax = dst_y_rsc, fill = "Density"), 
      color = col_dst
    ) +
    geom_hline(yintercept = 1, color = "gray") +
    geom_line(
      aes(x = platelet, y = lower_95, color = "95% CI"), 
      point, 
      linetype = 2
    ) +
    geom_line(
      aes(x = platelet, y = upper_95, color = "95% CI"), 
      point, 
      linetype = 2
    ) +
    geom_line(
      aes(x = platelet, y = HR, color = "Hazard Ratio"), 
      point, 
      linewidth = 0.9, 
      linetype = 1
    ) +
    annotate(
      "text", 
      x = max(point$platelet), 
      y = 0.3,
      # y = ceiling(max(dst_y_rsc)), 
      label = ifelse(
        p < 0.001, 
        "p-value < 0.001", 
        paste("p-value", "=", round(p, 3))
      ),
      hjust = 0.8
    ) +
    scale_x_continuous(
      # expand = c(0, 0),
      name = expression(paste("Platelet Count (", 10^9, "/L)"))
    ) +
    scale_y_continuous(
      # expand = c(0, 0), 
      breaks = seq(0, ceiling(hr_lmt[2]), by = 1),
      name = "Hazard Ratio", 
      sec.axis = sec_axis(
        trans = ~ rescale(., dst_lmt, c(0, max(.))), 
        name = "Density of Platelet Count"
      )
    ) +
    scale_color_manual(
      breaks = c("Hazard Ratio", "95% CI"), 
      values = c(col_hr, col_ci), 
      guide = guide_legend(
        title = NULL, 
        override.aes = list(linetype = c(1, 2), linewidth = 0.6)
      )
    ) +
    scale_fill_manual(
      name = NULL, 
      breaks = c("Density"),
      values = c(alpha(col_dst, 0.8))
    ) +
    coord_cartesian(
      ylim = c(0, ceiling(hr_lmt[2]) + 0.1)
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black")
    )
  
  return(p1)
}

# calculation -------------------------------------------------------------

vars <- c(
  "fu_time", "cancer_death",
  "platelet", "age",
  "sex",
  "asprin", "smoking_status", "alcohol_status", "bmi", "TDI", "race"
)

data <- extract_Cox_data(lagtime = 0, vars = vars)

hr_smth <- lapply(data, function(x) try(cal_hr(x), silent = TRUE))

# plotting ----------------------------------------------------------------

plot_list <- lapply(hr_smth, function(x) try(geom_hr(x), silent = TRUE)) %>% 
  extract(., which(sapply(., function(x) typeof(x) == "list")))

mplot <- plot_grid(
  plotlist = plot_list, 
  align = "none", 
  ncol = sqrt(length(plot_list)) %>% ceiling(), 
  labels = names(plot_list) %>% 
    gsub("_", " ", ., fixed = TRUE) %>% 
    gsub("E ", "Excluding ", ., fixed = TRUE),
  label_fontface = "plain",
  hjust = 0,
  vjust = 1.1,
  label_x = 0.1
)

save_plot("04/rcs_plot_lag0.pdf", mplot, base_height = 30)
