# env settings -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(scales)
library(patchwork)
library(ggtext)

load("02/density_plot_of_lag_time_in_pan-cancer.RData")
load("04/rcs_plot_pan_os.RData")
load("10/scatter_mr_pan_os.RData")
load("14/km_tpm4_os.RData")

dir.create("17", FALSE)

# plotting ---------------------------------------------------------------------

p <- wrap_elements(
  lag_plot_pan &
    theme(axis.text = element_text(color = "black"))
) +
  wrap_elements(
    p_pan_os +
      theme_classic() +
      theme(
        axis.text = element_text(color = "black"),
        plot.title = element_blank(),
        plot.caption = element_blank(),
        legend.position = "top"
      )
  ) +
  wrap_elements(
    scatter_pan_os +
      theme(
        axis.text = element_text(color = "black"),
        plot.title = element_blank(),
        legend.position = "top"
      )
  ) +
  wrap_elements(
    (
      tpm4_os[["plot"]] +
        theme_classic() +
        theme(
          axis.text = element_text(color = "black"),
          axis.title.x = element_blank()
        )
    ) + (
      tpm4_os[["table"]] +
        theme_classic() +
        theme(
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_markdown(color = c("#BC3C29FF", "#0072B5FF"))
        )
    ) +
      plot_layout(ncol = 1, heights = c(8, 2)) &
      theme(legend.position = "top")
  ) +
  plot_layout(ncol = 2, byrow = TRUE, widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

# plots saving -----------------------------------------------------------------

ggsave(
  "17/plots_combination.pdf",
  p,
  width = 13,
  height = 12
)
