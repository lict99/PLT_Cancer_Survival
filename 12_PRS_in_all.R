# env settings -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(patchwork)

load("01/UKB_all_info.RData")
load("11/PRS.RData")

source("functions/PRS_plot.R", local = TRUE)

dir.create("12", FALSE)

# calculation ------------------------------------------------------------------

## PRS for all participants in the UK biobank
prs_all <- merge(
  UKb_baseline[, c("eid", "platelet")],
  prs,
  by = "eid"
) %>%
  na.omit()

# plotting ---------------------------------------------------------------------

## weighted PRS plot
prs_w_plot <- geom_prs(
  data = prs_all,
  annotate = "All participants in UK Biobank",
  x_col = "prs_w"
)

## unweighted PRS plot
prs_u_plot <- geom_prs(
  data = prs_all,
  annotate = "All participants in UK Biobank",
  x_col = "prs_u"
)

# plots saving -----------------------------------------------------------------

ggsave(
  "12/weighted_PRS_in_all.pdf",
  prs_w_plot,
  width = 9.5,
  height = 9.5,
  units = "in"
)

ggsave(
  "12/unweighted_PRS_in_all.pdf",
  prs_u_plot,
  width = 9.5,
  height = 9.5,
  units = "in"
)
