# env settings ------------------------------------------------------------

library(magrittr)
library(survival)
library(ggplot2)
library(ggsci)
library(scales)
library(cowplot)

load("01/whole_cancer_data_for_Cox.RData")
load("09/genotype_with_cancer.RData")
load("09/pc_snp_summary.RData")
load("10/MR_harmonised_data.RData")
load("00/cancer_ICD_codes_with_attr.RData")
load("00/functions.RData")

dir.create("11", FALSE)

# PRS calculation ---------------------------------------------------------

## extract IV information
iv_info <- lapply(
  mr_data,
  function(x) {
    subset(x, mr_keep == TRUE) %>%
      merge(
        snp_smr[, c("affy", "snp")],
        by.x = "SNP",
        by.y = "snp",
        all.x = TRUE
      )
  }
)

## calculate SNP effect allele score
snp_score <- data.frame(eid = snp_ind_ca$eid) %>%
  (function(x) {
    for (i in iv_info$All_sites$SNP) {
      affy <- snp_smr[snp_smr$snp == i, "affy"]
      col_geno <- snp_ind_ca[, affy]
      snp_info <- subset(iv_info$All_sites, SNP == i)
      eff <- snp_info[, "effect_allele.exposure"]
      oth <- snp_info[, "other_allele.exposure"]
      homo_oth <- paste(oth, oth)
      homo_eff <- paste(eff, eff)
      heto <- c(paste(oth, eff), paste(eff, oth))
      for (j in seq_along(col_geno)) {
        if (!is.na(col_geno[j])) {
          if (col_geno[j] == homo_oth) {
            col_geno[j] <- 0
          } else if (col_geno[j] == heto[1]) {
            col_geno[j] <- 1
          } else if (col_geno[j] == heto[2]) {
            col_geno[j] <- 1
          } else if (col_geno[j] == homo_eff) {
            col_geno[j] <- 2
          } else {
            col_geno[j] <- NA
          }
        }
      }
      x[, i] <- as.numeric(col_geno)
    }
    return(x)
  })()

## un-weighted and weighted PRS for various cancers
prs <- data.frame(
  eid = snp_score$eid,
  prs_u = rowSums(snp_score[, -1]),
  prs_w = (function() {
    if (identical(names(snp_score[, -1]), iv_info$All_sites$SNP)) {
      a <- as.matrix(snp_score[, -1])
      b <- matrix(iv_info$All_sites$beta.exposure, ncol = 1)
    } else {
      stop("Not indentical in names!")
    }
    return(a %*% b)
  })()
)

prs_cancer <- extract_Cox_data(
  vars = c("eid", "fu_time", "cancer_death", "age", "sex", "platelet"),
  lagtime = 0
) %>%
  lapply(
    function(x) {
      merge(x, prs, by = "eid", all.x = TRUE) %>%
        na.omit()
    }
  )

# plotting ----------------------------------------------------------------

## function
gg_prs <- function(data, x_col, site) {
  # Calculation
  colors <- pal_nejm("default")(8) # %T>% show_col()
  col_d <- colors[5]
  col_l <- colors[1]
  mrg <- 2

  cor <- cor(data[, x_col], data$platelet, method = "pearson") %>%
    round(3)
  cor.p <- cor.test(
    data[, x_col],
    data$platelet,
    method = "pearson",
    alternative = "two.sided"
  ) %>%
    extract2("p.value") %>%
    round(3)

  if (x_col == "prs_w") {
    fit.smr <- lm(platelet ~ prs_w, data) %>%
      summary()
  } else if (x_col == "prs_u") {
    fit.smr <- lm(platelet ~ prs_u, data) %>%
      summary()
  }
  r.sq <- fit.smr$r.squared %>% round(3)
  f <- fit.smr$fstatistic[1] %>% round(3)
  fv1 <- fit.smr$fstatistic[2]
  fv2 <- fit.smr$fstatistic[3]
  fit.p <- pf(
    fit.smr$fstatistic[1],
    fit.smr$fstatistic[2],
    fit.smr$fstatistic[3],
    lower.tail = FALSE
  ) %>%
    round(3)

  # ggplot
  p1 <- ggplot() +
    geom_density(
      aes(x = data[, x_col]), 
      fill = alpha(col_d, 0.8), 
      color = col_d
    ) +
    labs(
      x = NULL,
      y = if (x_col == "prs_w") {
        "Density of Weighted PRS"
      } else if (x_col == "prs_u") {
        "Density of Unweighted PRS"
      } else {
        stop("Invalid label!")
      }
    ) +
    theme_classic() +
    theme(
      plot.margin = margin(mrg, 0, 0, mrg),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.length.x = unit(0, "pt")
    )

  p2 <- ggplot() +
    geom_blank() +
    theme_void() +
    theme(
      plot.margin = margin(mrg, mrg, 0, 0),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.ticks.length = unit(0, "pt")
    )

  p3 <- ggplot() +
    geom_point(
      aes(x = data[, x_col], y = data$platelet),
      alpha = 0.15,
      size = 0.1
    ) +
    geom_smooth(
      aes(x = data[, x_col], y = data$platelet),
      method = "lm",
      formula = "y ~ x",
      color = col_l
    ) +
    xlab(
      if (x_col == "prs_w") {
        "Weighted PRS"
      } else if (x_col == "prs_u") {
        "Unweighted PRS"
      } else {
        stop("Invalid label!")
      }
    ) +
    ylab(expression(paste("Platelet Count (", 10^9, "/L)", sep = ""))) +
    annotate(
      "text",
      x = max(data[, x_col]),
      y = max(data$platelet),
      label =
        paste(
          "cancer site: ", site,
          "\n",
          "sample size = ", nrow(data),
          "\n",
          "pearson r = ", cor, "; p-value = ", cor.p,
          "\n",
          "R-squared", " = ", r.sq,
          "\n",
          "F(", fv1, ", ", fv2, ") = ", f, "; p-value = ", fit.p,
          sep = ""
        ),
      hjust = 1,
      vjust = 1
    ) +
    theme_classic() +
    theme(
      plot.margin = margin(0, 0, mrg, mrg)
    )

  p4 <- ggplot(data, aes(platelet)) +
    geom_density(fill = alpha(col_d, 0.8), color = col_d) +
    labs(x = NULL, y = "Density of Platelet Count") +
    coord_flip() +
    theme_classic() +
    theme(
      plot.margin = margin(0, mrg, mrg, 0),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.length.y = unit(0, "pt")
    )

  plot_grid(
    p1, p2, p3, p4,
    nrow = 2,
    align = "hv",
    rel_heights = c(3, 10),
    rel_widths = c(10, 3)
  )
}

gg_prs_loop <- function(datalist, x_col) {
  plotlist <- list()
  for (i in names(datalist)) {
    plotlist[[i]] <- do.call(
      gg_prs,
      list(data = datalist[[i]], site = i, x_col = x_col)
    )
  }
  plotlist
}

## plot PRS
prs_w_plot <- gg_prs_loop(
  datalist = prs_cancer,
  x_col = "prs_w"
)

prs_u_plot <- gg_prs_loop(
  datalist = prs_cancer,
  x_col = "prs_u"
)

# PRS Cox -----------------------------------------------------------------

prs_w_cox <- lapply(
  prs_cancer,
  function(x) {
    if ("sex" %in% colnames(x)) {
      coxph(
        Surv(fu_time, cancer_death == 1) ~ prs_w + age + sex,
        data = x,
        singular.ok = TRUE,
        model = TRUE
      ) %>% 
        summary()
    } else {
      coxph(
        Surv(fu_time, cancer_death == 1) ~ prs_w + age,
        data = x,
        singular.ok = TRUE,
        model = TRUE
      ) %>% 
        summary()
    }
  }
)

prs_u_cox <- lapply(
  prs_cancer,
  function(x) {
    if ("sex" %in% colnames(x)) {
      coxph(
        Surv(fu_time, cancer_death == 1) ~ prs_u + age + sex,
        data = x,
        singular.ok = TRUE,
        model = TRUE
      ) %>% 
        summary()
    } else {
      coxph(
        Surv(fu_time, cancer_death == 1) ~ prs_u + age,
        data = x,
        singular.ok = TRUE,
        model = TRUE
      ) %>% 
        summary()
    }
  }
)

# data saving -------------------------------------------------------------

save.image(file = "11/image.RData")
save(iv_info, file = "11/IV_info.RData")
save(gg_prs, gg_prs_loop, file = "11/functions_gg_prs.RData")
