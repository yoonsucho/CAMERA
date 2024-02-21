#' @description
#' This method performs a IVW analysis on the harmonised data. It fits two linear models, one considering all populations and another considering each population separately. It then performs a fixed effects meta-analysis to estimate the heterogeneity across different populations. The results are stored in the `mrres` attribute of the `CAMERA` object.
#'
#' @param dat A data frame containing the harmonised data. It should have the columns `beta.y`, `beta.x`, `se.y`, and `pops`. If not provided, the method uses the `harmonised_dat` attribute of the `CAMERA` object.
#'
#' @return A list containing the results of the analysis. The list includes the coefficients of the fitted models, and the results of the heterogeneity analysis.
CAMERA$set("public", "cross_estimate", function(dat=self$harmonised_dat) {
  mod1 <- lm(beta.y ~ -1 + beta.x, data=dat, weight=1/dat$se.y^2)
  mod2 <- lm(beta.y ~ -1 + beta.x:as.factor(pops), data=dat, weight=1/dat$se.y^2)
  smod1 <- summary(mod1)
  smod2 <- summary(mod2)
  modcomp <- anova(mod1, mod2)

  res <- list()
  res$coefficients <- dplyr::bind_rows(
    smod1$coefficients %>% dplyr::as_tibble() %>% dplyr::mutate(pops="All"),
    smod2$coefficients %>% dplyr::as_tibble() %>% dplyr::mutate(pops=levels(as.factor(dat$pops)))
  ) %>%
    dplyr::select(pops, dplyr::everything())
  heterogeneity <- fixed_effects_meta_analysis(res$coefficients$Estimate[-1], res$coefficients$`Std. Error`[-1])
  heterogeneity$Qj <- dplyr::tibble(
    pops = res$coefficients$pops[-1],
    Qj = heterogeneity$Qj,
    Qjpval = heterogeneity$Qjpval
  )
  heterogeneity <- heterogeneity[names(heterogeneity) != "Qjpval"]
  res$coefficients <- dplyr::left_join(res$coefficients, heterogeneity$Qj, by="pops")
  res$coefficients$Qdf <- 1
  res$coefficients$Qj[res$coefficients$pops=="All"] <- heterogeneity$Q
  res$coefficients$Qjpval[res$coefficients$pops=="All"] <- heterogeneity$Qpval
  res$coefficients$Qdf[res$coefficients$pops=="All"] <- heterogeneity$Qdf
  self$mrres <- res$coefficients
  return(self$mrres)
})

#' @description
#' Plot the results from `cross_estimate`
#'
#' @param dat A data frame containing the harmonised data. It should have the columns `beta.y`, `beta.x`, `se.y`, and `pops`. If not provided, the method uses the `harmonised_dat` attribute of the `CAMERA` object.
#'
#' @return A list containing the results of the analysis. The list includes the coefficients of the fitted models, and the results of the heterogeneity analysis.
CAMERA$set("public", "plot_cross_estimate", function(res=self$mrres, qj_alpha=0.05) {
  est <- res$coefficients
  est$what <- "Pops"
  est$what[est$pops=="All"] <- "All"
  p <- ggplot2::ggplot(est, ggplot2::aes(x=Estimate, y=pops)) +
  ggplot2::geom_point(ggplot2::aes(colour=Qjpval < qj_alpha)) +
  ggplot2::geom_errorbarh(ggplot2::aes(colour=Qjpval < qj_alpha, xmin=Estimate - 1.96 * `Std. Error`, xmax=Estimate + 1.96 * `Std. Error`), height=0.1) +
  ggplot2::facet_grid(what ~ ., scale="free_y", space="free_y") +
  ggplot2::geom_vline(xintercept=0, linetype="dotted") +
  ggplot2::labs(y="", colour=paste0("Heterogeneity\npval < ", qj_alpha))
  return(p)
})

#' Identify blown up estimates
#' 
#' Sometimes estimates appear unstable. They are likely unreliable and best to not use for heterogeneity analyses etc. 
#' 
#' @param b Vector of betas
#' @param se Vector of SEs
#' @param infl Inflation factor - how much larger is the estimate than the estimate of the tightest SE
#' 
#' @export 
#' @return index of betas to remove
identify_blownup_estimates <- function(b, se, infl) {
  semin <- which.min(se)
  abs(b) > infl*abs(b[semin])
}

#' Perform fixed effects meta analysis for one association
#' 
#' @param beta_vec
#' @param se_vec
#' @param infl Inflation factor - how much larger is the estimate than the estimate of the tightest SE - for use in removing unreliable estimates
#' 
#' @export
#' @return list of results
fixed_effects_meta_analysis <- function(beta_vec, se_vec, infl=10000) {
    ind <- identify_blownup_estimates(beta_vec, se_vec, infl)
    beta_vec[ind] <- NA
    se_vec[ind] <- NA
    w <- 1 / se_vec^2
    beta <- sum(beta_vec * w, na.rm=T) / sum(w, na.rm=T)
    se <- sqrt(1 / sum(w, na.rm=T))
    pval <- pnorm(abs(beta / se), lower.tail = FALSE)
    Qj <- w * (beta-beta_vec)^2
    Q <- sum(Qj, na.rm=T)
    Qdf <- sum(!is.na(beta_vec))-1
    if(Qdf == 0) Q <- 0
    Qjpval <- pchisq(Qj, 1, lower.tail=FALSE)
    Qpval <- pchisq(Q, Qdf, lower.tail=FALSE)
    return(list(beta=beta, se=se, pval=pval, Q=Q, Qdf=Qdf, Qpval=Qpval, Qj=Qj, Qjpval=Qjpval))
}
