# env settings -----------------------------------------------------------------

library(magrittr)
library(table1)
library(openxlsx)
library(ggplot2)
library(patchwork)

load(file = "01/whole_cancer_data_for_Cox.RData")
load(file = "00/cancer_ICD_codes_with_attr.RData")
load(file = "00/cancer_names.RData")

source("functions/Cox_regression.R", local = TRUE)

dir.create("02", FALSE)

# statistics -------------------------------------------------------------------

## variables to be considered in this study
vars <- c(
  "platelet", "age", "sex", "bmi", "TDI", "aspirin",
  "smoking_status", "alcohol_status", "race",
  "fu_time", "fu_event", "lag_time"
)

## produce density plots of lag time
lag_plot_pan <- extract_Cox_data(
  data_list = whole_cancer_data[["OS"]],
  vars = vars,
  lagtime = c(0, Inf)
) %>%
  {
    data <- transform(
      extract2(., "All_sites"),
      lag_time = as.numeric(lag_time) %>% divide_by(365.25)
    )
    range_lag <- range(data[["lag_time"]])
    p1 <- ggplot(data = data, aes(x = lag_time)) +
      geom_histogram(
        aes(y = after_stat(density)),
        color = "#7876B1FF",
        fill = alpha("#7876B1FF", 0.8),
        binwidth = 0.5
      ) +
      geom_density(
        color = "#BC3C29FF"
      ) +
      labs(
        y = "Density of lag time"
      ) +
      coord_cartesian(
        xlim = c(floor(range_lag[[1]]), ceiling(range_lag[[2]]))
      ) +
      theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"),
        plot.margin = margin()
      )

    p2 <- ggplot(data = data, aes(y = lag_time)) +
      stat_boxplot(
        geom = "errorbar",
        width = 0.1,
        color = "#BC3C29FF"
      ) +
      geom_boxplot(
        width = 0.2,
        color = "#BC3C29FF"
      ) +
      labs(y = "Lag time (years)") +
      coord_flip(
        ylim = c(floor(range_lag[[1]]), ceiling(range_lag[[2]]))
      ) +
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        plot.margin = margin()
      )

    p <- p1 + p2 + plot_layout(ncol = 1, heights = c(9, 1))
    p
  }

## produce tables
tables <- lapply(
  whole_cancer_data,
  function(x) {
    extract_Cox_data(
      data_list = x,
      vars = vars,
      lagtime = c(-Inf, Inf)
    ) %>%
      lapply(
        function(y) {
          transform(
            y,
            aspirin = factor(aspirin, levels = c("YES", "NO")),
            race = factor(race, levels = c("British", "Others")),
            lag_group = ifelse(
              lag_time >= 0,
              "lag_time >= 0",
              "lag_time < 0"
            ) %>%
              factor(),
            lag_time = as.numeric(lag_time) %>% divide_by(365.25),
            fu_time = as.numeric(fu_time) %>% divide_by(365.25),
            fu_event = factor(
              fu_event,
              levels = c(1, 0),
              labels = c("Death", "Censoring")
            )
          ) %>%
            table1(
              ~ . | lag_group,
              data = .,
              render.continuous = function(x) {
                with(
                  stats.default(x),
                  c(
                    "",
                    `Median (IQR)` = sprintf(
                      "%.2f (%.2f-%.2f)", MEDIAN, Q1, Q3
                    ),
                    `Min - Max` = sprintf(
                      "%.2f-%.2f", MIN, MAX
                    )
                  )
                )
              },
              overall = FALSE,
              extra.col = NULL
            ) %>%
            as.data.frame()
        }
      )
  }
)

## calculate number and proportion
## all participants with cancer including eligible and ineligible patients
nprop <- list(
  lag_no_limit = extract_Cox_data(
    data_list = whole_cancer_data[["OS"]],
    vars = c("eid"),
    lagtime = c(-Inf, Inf)
  ) %>%
    {
      n <- sapply(., nrow)
      n_max <- max(n)
      prop <- sapply(., function(x) nrow(x) / n_max)
      data.frame(
        cancer = unlist(cancer_names[names(.)]),
        N = n,
        proportion = paste0(sprintf("%.1f", prop * 100), "%")
      )
    },
  lag_0 = extract_Cox_data(
    data_list = whole_cancer_data[["OS"]],
    vars = c("eid"),
    lagtime = c(0, Inf)
  ) %>%
    {
      n <- sapply(., nrow)
      n_max <- max(n)
      prop <- sapply(., function(x) nrow(x) / n_max)
      data.frame(
        cancer = unlist(cancer_names[names(.)]),
        N = n,
        proportion = paste0(sprintf("%.1f", prop * 100), "%")
      )
    }
)

# results saving ---------------------------------------------------------------

write.xlsx(tables[["OS"]], file = "02/table1_for_OS.xlsx")
write.xlsx(tables[["CSS"]], file = "02/table1_for_CSS.xlsx")
write.xlsx(nprop, file = "02/number_and_proportion.xlsx")

# data saving ------------------------------------------------------------------

save(
  lag_plot_pan,
  file = "02/density_plot_of_lag_time_in_pan-cancer.RData",
  compress = FALSE
)
