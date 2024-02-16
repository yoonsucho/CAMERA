#' Estimate expected vs observed replication of effects between discovery and replication datasets
#' 
#' Taken from Okbay et al 2016. Under the assumption that all discovery effects are unbiased, what fraction of associations would replicate in the replication dataset, given the differential power of the discovery and replication datasets.
#' Uses standard error of the replication dataset to account for differences in sample size and distribution of independent variable
#' 
#' @param b_disc Vector of discovery betas
#' @param b_rep Vector of replication betas
#' @param se_disc Vector of discovery standard errors
#' @param se_rep Vector of replication standard errors
#' @param alpha Nominal replication significance threshold
#' 
#' @return List of results
#' - res: aggregate expected replication rate vs observed replication rate
#' - variants: per variant expected replication rates
prop_overlap <- function(b_disc, b_rep, se_disc, se_rep, alpha) {
  p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))
  p_rep <- pnorm(abs(b_rep) / se_rep, lower.tail = FALSE)
  res <- tibble::tibble(
    nsnp = length(b_disc),
    metric = c("Sign", "Sign", "P-value", "P-value"),
    datum = c("Expected", "Observed", "Expected", "Observed"),
    value = c(sum(p_sign, na.rm = TRUE), sum(sign(b_disc) == sign(b_rep)), sum(p_sig, na.rm = TRUE), sum(p_rep < alpha, na.rm = TRUE))
  ) %>%
    dplyr::group_by(metric) %>%
      dplyr::do({
        x <- .
        if(.$nsnp[1] > 0) {
          bt <- binom.test(
            x=.$value[.$datum == "Observed"], 
            n=.$nsnp[1], 
            p=.$value[.$datum == "Expected"] / .$nsnp[1]
          )$p.value
          x$pdiff <- bt
        }
        x
      })
  return(list(res = res, variants = dplyr::tibble(sig = p_sig, sign = p_sign, )))
}

#' @description
#' The function evaluates heterogeneity in the association of selected instruments and the exposure/outcome between the populations. Heterogeneity (Q statistics) is calcuated based on an IVW or simple MODE MR estimator. The instruments can be identified using "Raw", "MaxZ", or fine-mapping (Susie, PAINTOR) methods.
#' @param instrument Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param alpha Statistical threshold to determine significance. Default is "bonferroni", which is eqaul to 0.05/number of the instruments.
#' @param method IVW or Simple MODE
#' @return Table of the result
CAMERA$set("public", "instrument_heterogeneity", function(instrument = self$instrument_raw, alpha = "bonferroni", method = "ivw", outlier_removal = FALSE) {
  if (alpha == "bonferroni") {
    alpha <- 0.05 / (nrow(instrument))
  }


  if (outlier_removal == TRUE) {
    d <- dplyr::inner_join(
      subset(instrument, id == self$exposure_ids[[1]]),
      subset(instrument, id == self$exposure_ids[[2]]),
      by = "rsid"
    )

    rd <- RadialMR::format_radial(d$beta.x, d$beta.y, d$se.x, d$se.y, d$rsid)
    outliers <- RadialMR::ivw_radial(rd, alpha, 3, 0.0001, FALSE)$outliers

    if (outliers[1] != "No significant outliers") {
      outliers <- outliers$SNP
      instrument <- subset(instrument, !rsid %in% outliers)
    } else {
      print("No significant outliers")
    }
  }

  if (method == "ivw") {
    if (!any(names(instrument) %in% c("beta.outcome"))) {
      o <- lapply(self$exposure_ids, function(i) {
        m <- subset(instrument, id == i & p < alpha)
        other_ids <- self$exposure_ids[!self$exposure_ids %in% i]

        o <- lapply(other_ids, function(j) {
          n <- subset(instrument, id == j & rsid %in% m$rsid)
          dat <- dplyr::inner_join(m, n, by = "rsid") %>%
            dplyr::select(SNP = rsid, x = beta.x, y = beta.y, xse = se.x, yse = se.y, xp = p.x, yp = p.y)
          res <- suppressMessages(TwoSampleMR::mr_ivw(dat$x, dat$y, dat$xse, dat$yse)) %>%
            {
              tibble::tibble(Reference = i, Replication = j, nsnp = length(unique(dat$SNP)), agreement = .$b, se = .$se, pval = .$pval, Q = .$Q, Q_pval = .$Q_pval, I2 = ((.$Q - length(dat$SNP)) / .$Q))
            }
          return(res)
        })
      })
    }

    if (any(names(instrument) %in% c("beta.outcome"))) {
      o <- lapply(self$outcome_ids, function(i) {
        m <- subset(instrument, id.outcome == i & pval.outcome < alpha)
        other_ids <- self$outcome_ids[!self$outcome_ids %in% i]

        o <- lapply(other_ids, function(j) {
          n <- subset(instrument, id.outcome == j & SNP %in% m$SNP)
          dat <- dplyr::inner_join(m, n, by = "SNP") %>%
            dplyr::select(SNP = SNP, x = beta.outcome.x, y = beta.outcome.y, xse = se.outcome.x, yse = se.outcome.y, xp = pval.outcome.x, yp = pval.outcome.y)
          res <- suppressMessages(TwoSampleMR::mr_ivw(dat$x, dat$y, dat$xse, dat$yse)) %>%
            {
              tibble::tibble(Reference = i, Replication = j, nsnp = length(unique(dat$SNP)), agreement = .$b, se = .$se, pval = .$pval, Q = .$Q, Q_pval = .$Q_pval, I2 = ((.$Q - length(dat$SNP)) / .$Q))
            }
          return(res)
        })
      })
    }
  }

  if (method == "simple_mode") {
    if (!any(names(instrument) %in% c("beta.outcome"))) {
      o <- lapply(self$exposure_ids, function(i) {
        m <- subset(instrument, id == i & p < alpha)
        other_ids <- self$exposure_ids[!self$exposure_ids %in% i]

        o <- lapply(other_ids, function(j) {
          n <- subset(instrument, id == j & rsid %in% m$rsid)
          dat <- dplyr::inner_join(m, n, by = "rsid") %>%
            dplyr::select(SNP = rsid, x = beta.x, y = beta.y, xse = se.x, yse = se.y, xp = p.x, yp = p.y)
          res <- suppressMessages(TwoSampleMR::mr_simple_mode(dat$x, dat$y, dat$xse, dat$yse)) %>%
            {
              tibble::tibble(Reference = i, Replication = j, nsnp = length(unique(dat$SNP)), agreement = .$b, se = .$se, pval = .$pval)
            }
          return(res)
        })
      })
    }

    if (any(names(instrument) %in% c("beta.outcome"))) {
      o <- lapply(self$outcome_ids, function(i) {
        m <- subset(instrument, id.outcome == i & pval.outcome < alpha)
        other_ids <- self$outcome_ids[!self$outcome_ids %in% i]

        o <- lapply(other_ids, function(j) {
          n <- subset(instrument, id.outcome == j & SNP %in% m$SNP)
          dat <- dplyr::inner_join(m, n, by = "SNP") %>%
            dplyr::select(SNP = SNP, x = beta.outcome.x, y = beta.outcome.y, xse = se.outcome.x, yse = se.outcome.y, xp = pval.outcome.x, yp = pval.outcome.y)
          res <- suppressMessages(TwoSampleMR::mr_simple_mode(dat$x, dat$y, dat$xse, dat$yse)) %>%
            {
              tibble::tibble(Reference = i, Replication = j, nsnp = length(unique(dat$SNP)), agreement = .$b, se = .$se, pval = .$pval)
            }
          return(res)
        })
      })
    }
  }
  print(o %>% dplyr::bind_rows())
})


# Once a set of instruments is chosen we can ask
# - what fraction of those primarily identified in pop1 replicate in pop2?
# - what fraction of those primarily identified in pop1 have the same sign as in pop2?
# - compare these to what is expected by chance under the hypothesis that the effect estimates are the same
# this will provide some evidence for whether lack of replication is due to (e.g.) GxE
#' @description
#' The function estimates what fraction of the instuments is expected to be replicated across the populations under the hypothesis that the effect estimates are the same.
#' @param instrument Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param alpha Statistical threshold to determine significance. Default is "bonferroni", which is eqaul to 0.05/number of the instruments.
#' @param winnerscurse Use this option to correct winners' curse bias.
#' @return Table of the results. Summary of the results available in x$instrument_specificity_summary.
CAMERA$set("public", "estimate_instrument_specificity", function(instrument, alpha = "bonferroni", winnerscurse = FALSE) {
  if (alpha == "bonferroni") {
    alpha <- 0.05 / nrow(instrument)
  }
  # for every exposure that has instruments
  # get its instruments
  # check for replication in every other exposure
  o <- lapply(self$exposure_ids, function(i) {
    m <- subset(instrument, id == i & p < alpha)
    other_ids <- self$exposure_ids[!self$exposure_ids %in% i]

    if (winnerscurse == TRUE) {
      wcm <- m %>% dplyr::select(rsid, beta, se)
      m <- winnerscurse::cl_interval(summary_data = wcm, alpha = alpha, conf_level = 0.95) %>%
        dplyr::mutate(se = (.$upper - .$lower) / 3.92) %>%
        dplyr::arrange(rsido) %>%
        dplyr::rename(beta = beta.cl3)
    }
    o <- lapply(other_ids, function(j) {
      message("Checking ", i, " against ", j)
      other <- subset(instrument, id == j)
      temp <- dplyr::inner_join(m, other, by = "rsid")
      o <- other %>%
        {
          prop_overlap(temp$beta.x, temp$beta.y, temp$se.x, temp$se.y, alpha)
        }
      o$res <- o$res %>%
        dplyr::mutate(discovery = i, replication = j) %>%
        dplyr::select(discovery, replication, dplyr::everything())
      o$variants <- o$variants %>%
        dplyr::mutate(
          discovery = i,
          replication = j,
          rsid = temp$rsid,
          p_rep = temp$p.y,
          sign_agreement = sign(temp$beta.x) == sign(temp$beta.y),
          distinct = (sig > 0.8 & p_rep > 0.1) | (sign > 0.8 & !sign_agreement)
        ) %>%
        dplyr::select(discovery, replication, dplyr::everything())
      return(o)
    })
    overall <- lapply(o, function(x) {
      x$res
    }) %>% dplyr::bind_rows()
    pervariant <- lapply(o, function(x) {
      x$variants
    }) %>% dplyr::bind_rows()
    return(list(overall = overall, pervariant = pervariant))
  })
  self$instrument_specificity_summary <- lapply(o, function(x) x$overall) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(nsnp > 0)
  self$instrument_specificity <- lapply(o, function(x) x$pervariant) %>% dplyr::bind_rows()
  return(self$instrument_specificity_summary %>% as.data.frame())
})

#' @description
#' The fuction explains what contributes to the replication of gene-trait association between the populations, considering LD structure.
#' @param instrument Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param ld LD matrix obtained by using \code{x$regional_ld_matrices()}
#' @return Table of the regression result
CAMERA$set("public", "replication_evaluation", function(instrument = self$instrument_raw, ld = self$ld_matrices) {
  instrument_pop <- instrument %>% dplyr::group_split(id)

  rep <- instrument_pop[[1]] %>% dplyr::select(rsid)
  rep$sign <- sign(instrument_pop[[1]]$beta) == sign(instrument_pop[[2]]$beta)
  rep$sig <- instrument_pop[[1]]$p < 5e-8 & instrument_pop[[2]]$p < 5e-8

  fam1 <- read.table(paste0(self$bfiles[[1]], ".fam"), header = FALSE)
  n1 <- length(unique(fam1[[1]]))

  ldsc_pop1 <- lapply(1:length(ld), function(i) {
    ld <- ld[[i]][[1]]
    r2 <- ((n1 - 1) / (n1 - 2) * (ld^2)) - (1 / (n1 - 2))
    l <- mean(r2[lower.tri(r2)]^2, diag = FALSE)
    return(l)
  }) %>% unlist()

  fam2 <- read.table(paste0(self$bfiles[[2]], ".fam"), header = FALSE)
  n2 <- length(unique(fam2[[1]]))

  ldsc_pop2 <- lapply(1:length(ld), function(i) {
    ld <- ld[[i]][[2]]
    r2 <- ((n2 - 1) / (n2 - 2) * (ld^2)) - (1 / (n2 - 2))
    l <- mean(r2[lower.tri(r2)]^2, diag = FALSE)
    return(l)
  }) %>% unlist()

  rep$delta_ld <- ldsc_pop1 - ldsc_pop2

  instrument_pop <- lapply(1:length(self$exposure_ids), function(i) {
    instrument_pop[[i]] %>% dplyr::mutate(maf = dplyr::if_else(.$eaf > 0.5, (1 - .$eaf), .$eaf, NA_real_))
  })

  maf_pop1 <- instrument_pop[[1]]$maf * (1 - instrument_pop[[1]]$maf)
  maf_pop2 <- instrument_pop[[2]]$maf * (1 - instrument_pop[[2]]$maf)

  rep$delta_maf <- maf_pop1 - maf_pop2

  res <- list()
  res[[1]] <- summary(glm(sign ~ delta_ld, data = rep, family = binomial(link = "logit")))$coefficients
  res[[2]] <- summary(glm(sig ~ delta_maf, data = rep, family = binomial(link = "logit")))$coefficients
  res[[3]] <- summary(glm(sign ~ delta_ld + delta_maf, data = rep, family = binomial(link = "logit")))$coefficients
  res[[4]] <- summary(glm(sig ~ delta_ld + delta_maf, data = rep, family = binomial(link = "logit")))$coefficients
  names(res)[1] <- c("replicated_sign_model1")
  names(res)[2] <- c("replicated_sig_model1")
  names(res)[3] <- c("replicated_sign_model2")
  names(res)[4] <- c("replicated_sig_model2")
  return(list(res))
})
