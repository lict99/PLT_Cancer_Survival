# the nonlinear MR plot ----
# modified from nlmr::fracpoly_figure()

fracpoly_figure2 <- function(beta, cov, x.min, x.max, family = "gaussian", d = 1, p_ML = NULL, p1_ML = NULL, p2_ML = NULL, ci = "model_se", frac_coef_boot = NULL, ref = mean(x), pref_x = "x", pref_x_ref = "x", pref_y = "y", ci_type = "overall", ci_quantiles = 10, breaks = NULL, p = NULL) {
  if (ci_type == "overall") {
    plot.data <- data.frame(x = runif(10000, x.min, x.max))
    plot.data.ref <- data.frame(x = ref, y = 0)
    if (d == 1) {
      if (p_ML == -1) {
        plot.data$yest <- beta * log(plot.data$x) - (beta * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((log(plot.data$x) - log(ref))^2 * cov[1, 1])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot) - reprow(log(ref) %*% t(frac_coef_boot), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p_ML != -1) {
        plot.data$yest <- beta * plot.data$x^(p_ML + 1) - beta * ref^(p_ML + 1)
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((plot.data$x^(p_ML + 1) - ref^(p_ML + 1))^2 * cov[1, 1])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p_ML + 1) %*% t(frac_coef_boot[, 1]) - reprow(ref^(p_ML + 1) %*% t(frac_coef_boot[, 1]), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
    }
    if (d == 2) {
      if (p1_ML == -1 & p2_ML == -1) {
        plot.data$yest <- beta[1] * log(plot.data$x) + beta[2] * log(plot.data$x) * log(plot.data$x) - (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((log(plot.data$x) - log(ref))^2 * cov[1, 1] + 2 * (log(plot.data$x) - log(ref)) * (log(plot.data$x) * log(plot.data$x) - log(ref) * log(ref)) * cov[1, 2] + (log(plot.data$x) * log(plot.data$x) - log(ref) * log(ref))^2 * cov[2, 2])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot[, 1]) + log(plot.data$x) * log(plot.data$x) %*% t(frac_coef_boot[, 2]) - reprow(log(ref) %*% t(frac_coef_boot[, 1]) + log(ref) * log(ref) %*% t(frac_coef_boot[, 2]), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML == -1 & p2_ML != -1 & p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * log(plot.data$x) + beta[2] * plot.data$x^(p2_ML + 1) - (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((log(plot.data$x) - log(ref))^2 * cov[1, 1] + 2 * (log(plot.data$x) - log(ref)) * (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1)) * cov[1, 2] + (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1))^2 * cov[2, 2])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- log(plot.data$x) %*% t(frac_coef_boot[, 1]) + plot.data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) - reprow(log(ref) %*% t(frac_coef_boot[, 1]) + ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 & p2_ML == -1 & p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) + beta[2] * log(plot.data$x) - (beta[1] * ref^(p1_ML + 1) + beta[2] * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] + 2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (log(plot.data) - log(ref)) * cov[1, 2] + (log(plot.data$x) - log(ref))^2 * cov[2, 2])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) + log(plot.data$x) %*% t(frac_coef_boot[, 2]) - reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) + log(ref) %*% t(frac_coef_boot[, 2]), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 & p2_ML != -1 & p1_ML == p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) + beta[2] * plot.data$x^(p2_ML + 1) * log(plot.data$x) - (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1) * log(ref))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] + 2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (plot.data$x^(p2_ML + 1) * log(plot.data$x) - ref^(p2_ML + 1) * log(ref)) * cov[1, 2] + (plot.data$x^(p2_ML + 1) * log(plot.data$x) - ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) + (plot.data$x^(p2_ML + 1) * log(plot.data$x)) %*% t(frac_coef_boot[, 2]) - reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) + (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
      if (p1_ML != -1 & p2_ML != -1 & p1_ML != p2_ML) {
        plot.data$yest <- beta[1] * plot.data$x^(p1_ML + 1) + beta[2] * plot.data$x^(p2_ML + 1) - (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
        if (ci != "bootstrap_per") {
          plot.data$yse <- sqrt((plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1))^2 * cov[1, 1] + 2 * (plot.data$x^(p1_ML + 1) - ref^(p1_ML + 1)) * (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1)) * cov[1, 2] + (plot.data$x^(p2_ML + 1) - ref^(p2_ML + 1))^2 * cov[2, 2])
          plot.data$lci <- plot.data$yest - 1.96 * plot.data$yse
          plot.data$uci <- plot.data$yest + 1.96 * plot.data$yse
        } else {
          boot <- plot.data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) + plot.data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) - reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) + ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]), n = nrow(plot.data))
          plot.data$lci <- rowQuantiles(boot, probs = 0.025)
          plot.data$uci <- rowQuantiles(boot, probs = 0.975)
        }
      }
    }
    if (family != "binomial") {
      figure <- ggplot(plot.data, aes(x = x))
      figure <- figure + geom_hline(aes(yintercept = 0), colour = "grey") + geom_line(aes(y = yest), color = "black") + geom_line(aes(y = lci), color = "grey") + geom_line(aes(y = uci), color = "grey") + theme_bw() + labs(x = pref_x, y = bquote(.(pref_y) ~ " [" ~ .(pref_x_ref)["ref"] ~ "=" ~ .(round(ref, 2)) ~ "]")) + theme(axis.title.x = element_text(vjust = 0.5, size = 20), axis.title.y = element_text(vjust = 0.5, size = 20), axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) + geom_point(aes(x = x, y = y), data = plot.data.ref, colour = "red", size = 4) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      if (!is.null(breaks)) {
        suppressMessages(figure <- figure + scale_y_continuous(breaks = breaks))
      }
    }
    # where I modified
    if (family == "binomial") {
      plot.data$yest <- exp(plot.data$yest)
      plot.data$uci <- exp(plot.data$uci)
      plot.data$lci <- exp(plot.data$lci)
      plot.data.ref$y <- exp(0)
      figure <- ggplot(plot.data, aes(x = x)) +
        geom_hline(yintercept = 1, colour = "grey") +
        geom_line(aes(y = lci, color = "ci"), linetype = 1) +
        geom_line(aes(y = uci, color = "ci"), linetype = 1) +
        geom_line(aes(y = yest, color = "y"), linewidth = 1) +
        geom_point(
          aes(x = x, y = y),
          data = plot.data.ref,
          colour = "#7876B1FF",
          size = 2
        ) +
        labs(
          x = expression("Platelet counts (" * 10^9 * "/L)"),
          y = "Odds ratio"
        ) +
        scale_y_continuous(
          trans = "log",
          breaks = c(0, 1, pretty(plot.data$uci))
        ) +
        scale_color_manual(
          breaks = c("y", "ci"),
          labels = c("Odds ratio", "95% confidence interval"),
          values = c("#BC3C29FF", "#6F99ADFF"),
          guide = guide_legend(
            title = NULL,
            override.aes = list(linetype = c(1, 1), linewidth = 0.4)
          )
        ) +
        annotate(
          "label",
          x = max(plot.data$x),
          y = 1 + 0.5,
          label = paste0(
            "italic('p')['non-linearity']=='",
            sprintf("%.3f", p),
            "'"
          ),
          parse = TRUE,
          hjust = 1,
          vjust = 0,
          label.size = 0
        ) +
        theme_classic() +
        theme(
          legend.position.inside = c(1, 1),
          legend.justification = c(1, 1),
          legend.direction = "horizontal"
        )
    }
  } else {
    # delete other ci_type
    stop("You need remove this function in `.GlobalEnv`!")
  }
  figure
}

# nonlinear MR result ----
# modified from nlmr::fracpoly_mr()

fracpoly_mr2 <- function(y, x, g, covar = NULL, family = "gaussian", q = 10, xpos = "mean", method = "FE", d = 1, pd = 0.05, ci = "model_se", nboot = 100, fig = F, ref = mean(x), pref_x = "x", pref_x_ref = "x", pref_y = "y", ci_type = "overall", ci_quantiles = 10, breaks = NULL) {
  ##### Error messages #####
  if (!(is.vector(y) & is.vector(x) & is.vector(g))) stop("either the outcome, exposure or instrument is not a vector")
  if (!is.null(covar)) {
    if (!is.data.frame(covar)) stop("covar has to be a data.frame")
  }
  if (!((is.numeric(y) | is.integer(y)) & (is.numeric(x) | is.integer(x)) & (is.numeric(g) | is.integer(g)))) stop("either the outcome, exposure or instrument is not numeric")
  if (any(x <= 1)) stop("fractional polynomial models require the exposure to be >>1")
  if (length(y) <= 1) stop("the outcome is less than or equal to a single value")
  if (!(length(y) == length(x) & length(y) == length(g)) | (if (!is.null(covar)) {
    (nrow(covar) != length(y))
  } else {
    FALSE
  })) {
    stop("the number of observations for the outcome, exposure, instrument and covariates are not all the same")
  }
  if (any(is.na(y)) | any(is.na(x)) | any(is.na(g)) | (if (!is.null(covar)) {
    any(is.na(covar))
  } else {
    FALSE
  })) {
    stop("there are missing values in either the outcome, exposure, instrument or covariates")
  }
  if (!(family == "gaussian" | family == "binomial")) stop('family has to be equal to either "gaussian" or "binomial"')
  if (family == "binomial") {
    if (any(!(y == 1 | y == 0))) stop('y has to be 0 or 1 if family is equal to "binomial"')
  }
  if ((length(y) / 10) < q) stop("the quantiles should contain at least 10 observations")
  if (!(xpos == "mean" | (xpos > 0 & xpos < 1))) stop("the position used to relate x to the localised average causal effect")
  if (!(d == 1 | d == 2 | d == "both")) stop('the degree has to be equal to 1, 2 or "both"')
  if (!(ci == "model_se" | ci == "bootstrap_se" | ci == "bootstrap_per")) stop('the confidence intervals must be one of "model_se", "bootstrap_se" and "bootstrap_per"')

  ##### Covariates #####
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, -1, drop = F]
    if (any(is.na(covar))) stop("there are missing values in the covariates")
  }

  ##### x0 (IV-Free) #####
  ivf <- iv_free(y = y, x = x, g = g, covar = covar, q = q, family = family)
  x0 <- ivf$x0
  xcoef <- ivf$xcoef
  x0q <- ivf$x0q

  ##### LACE #####
  loc <- lace(y = y, x = x, g = g, covar = covar, q = q, x0q = x0q, xc_sub = TRUE, family = family, xpos = xpos)
  coef <- loc$coef / xcoef
  coef_se <- loc$coef_se / xcoef
  xmean <- loc$xmean
  xcoef_sub <- loc$xcoef_sub
  xcoef_sub_se <- loc$xcoef_sub_se

  ##### Test of IV-exposure assumption #####
  p_het <- 1 - pchisq(rma(xcoef_sub, vi = (xcoef_sub_se)^2)$QE, df = (q - 1))
  p_het_trend <- rma.uni(xcoef_sub ~ xmean, vi = xcoef_sub_se^2, method = method)$pval[2]

  ##### Best-fitting fractional polynomial #####
  fracpb <- fracpoly_best(coef = coef, coef_se = coef_se, xmean = xmean, d = d, pd = pd, method = method)
  model <- fracpb$model
  p_ML <- fracpb$p_ML
  p1_ML <- fracpb$p1_ML
  p2_ML <- fracpb$p2_ML
  fp_p <- fracpb$fp_p
  fp_d12_p <- fracpb$fp_d12_p
  d <- fracpb$d

  ##### Non-linearity tests #####
  p_quadratic <- rma(coef ~ xmean, (coef_se)^2, method = "FE")$pval[2]
  p_Q <- 1 - pchisq(rma(coef, vi = (coef_se)^2)$QE, df = (q - 1))

  ##### Bootstrap #####
  if (ci == "bootstrap_per" | ci == "bootstrap_se") {
    frac_coef_boot <- fracpoly_boot(y = y, x = x, g = g, covar = covar, q = q, x0q = x0q, xcoef = xcoef, family = family, xpos = xpos, method = method, nboot = nboot, d = d, p_ML = p_ML, p1_ML = p1_ML, p2_ML = p2_ML)
  } else {
    frac_coef_boot <- NULL
  }

  ##### Results #####
  if (d == 1) {
    powers <- p_ML + 1
  } else {
    powers <- c((p1_ML + 1), (p2_ML + 1))
  }
  beta <- as.numeric(model$b)
  if (ci == "model_se") {
    cov <- model$vb
    se <- model$se
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * pnorm(-abs(beta / se))
  }
  if (ci == "bootstrap_se") {
    cov <- var(frac_coef_boot)
    se <- sqrt(diag(cov))
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * pnorm(-abs(beta / se))
  }
  if (ci == "bootstrap_per") {
    if (d == 1) {
      se <- NA
      lci <- quantile(frac_coef_boot, probs = 0.025)
      uci <- quantile(frac_coef_boot, probs = 0.975)
      pval <- NA
    }
    if (d == 2) {
      se <- rep(NA, 2)
      lci <- NULL
      uci <- NULL
      pval <- NULL
      lci[1] <- quantile(frac_coef_boot[, 1], probs = 0.025)
      lci[2] <- quantile(frac_coef_boot[, 2], probs = 0.025)
      uci[1] <- quantile(frac_coef_boot[, 1], probs = 0.975)
      uci[2] <- quantile(frac_coef_boot[, 2], probs = 0.975)
      pval <- rep(NA, 2)
    }
  }
  lci <- as.numeric(lci)
  uci <- as.numeric(uci)
  if (ci == "model_se") {
    nboot <- NA
  }

  ##### Figure #####
  # where I modified
  if (fig == T) {
    if (d == 1) {
      figure <- fracpoly_figure2(beta = beta, cov = cov, x.min = min(x), x.max = max(x), family = family, d = d, p_ML = p_ML, ci = ci, frac_coef_boot = frac_coef_boot, ref = ref, pref_x = pref_x, pref_x_ref = pref_x_ref, pref_y = pref_y, ci_type = ci_type, ci_quantile = ci_quantile, breaks = breaks, p = fp_p)
    }
    if (d == 2) {
      figure <- fracpoly_figure2(beta = beta, cov = cov, x.min = min(x), x.max = max(x), family = family, d = d, p1_ML = p1_ML, p2_ML = p2_ML, ci = ci, frac_coef_boot = frac_coef_boot, ref = ref, pref_x = pref_x, pref_x_ref = pref_x_ref, pref_y = pref_y, ci_type = ci_type, ci_quantile = ci_quantile, breaks = breaks, p = fp_p)
    }
  }

  ##### Return #####
  model <- as.matrix(data.frame(q = q, xpos = xpos, ci_type = ci, nboot = nboot))
  coefficients <- as.matrix(data.frame(beta = beta, se = se, lci = lci, uci = uci, pval = pval))
  rownames(coefficients) <- powers
  if (nrow(coefficients) == 2) {
    if (powers[1] == powers[2]) {
      rownames(coefficients) <- c(powers[1], paste0("log ", powers[2]))
    }
  }
  loc <- as.matrix(data.frame(beta = (coef), se = (abs(coef_se)), lci = (coef - 1.96 * (abs(coef_se))), uci = (coef + 1.96 * (abs(coef_se))), pval = (2 * pnorm(-abs(coef / coef_se)))))
  rownames(loc) <- 1:nrow(loc)
  xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
  rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
  p_tests <- as.matrix(data.frame(fp_d1_d2 = fp_d12_p, fp = fp_p, quad = p_quadratic, Q = p_Q))
  p_heterogeneity <- as.matrix(data.frame(Q = p_het, trend = p_het_trend))
  if (fig == F) {
    results <- list(n = length(y), model = model, powers = powers, coefficients = coefficients, lace = loc, xcoef = xcoef_quant, p_tests = p_tests, p_heterogeneity = p_heterogeneity)
  } else {
    results <- list(n = length(y), model = model, powers = powers, coefficients = coefficients, lace = loc, xcoef = xcoef_quant, p_tests = p_tests, p_heterogeneity = p_heterogeneity, figure = figure)
  }
  class(results) <- "fracpoly_mr"
  return(results)
}
