# env settings ----

library(magrittr)
library(ggplot2)
library(patchwork)

load("00/cancer_names.RData")
load("03/platelet_Cox.RData")

dir.create("07", FALSE)

# data ----
