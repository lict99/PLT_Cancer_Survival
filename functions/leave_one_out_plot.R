# Leave-one-out sensitivity analysis plot created by ggplot2

geom_loo <- function(loo) {
  if (sum(!grepl("All", loo$SNP)) < 3) {
    stop("Number of SNPs < 3!")
  }
  loo$or <- exp(loo$b)
  loo$up <- exp(loo$b + 1.96 * loo$se)
  loo$low <- exp(loo$b - 1.96 * loo$se)
  loo <- rbind(
    subset(loo, SNP != "All") %>%
      extract(order(use_series(., or)), ),
    subset(loo, SNP == "All")
  )
  loo$index <- seq(nrow(loo) + 1, 1)[-(nrow(loo))]
  p <- ggplot(loo, aes(y = index, x = or)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_hline(yintercept = 2, color = "gray") +
    geom_pointrange(aes(xmin = low, xmax = up)) +
    geom_point(
      aes(x = or[length(or)], y = 1),
      color = "#BC3C29FF",
      size = 2.5
    ) +
    scale_y_continuous(breaks = loo$index, labels = loo$SNP) +
    labs(
      y = "",
      x = "Odds ratio (95% confidence interval)",
      title = loo$id.outcome,
      caption = paste(
        "Leave-one-out sensitivity analysis",
        "using inverse variance weighted method",
        "\n",
        "for", loo$exposure, "on",
        switch(unique(loo$outcome),
          "OS" = "overall survival",
          "CSS" = "cancer-specific survival",
          stop()
        )
      )
    ) +
    coord_cartesian(ylim = c(0, nrow(loo))) +
    theme_classic() +
    theme(
      plot.caption = element_text(hjust = 0.5)
    )
  p
}
