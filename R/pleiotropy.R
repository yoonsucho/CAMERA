CAMERA$set("public", "pleiotropy", function(harmonised_dat = self$harmonised_dat) {
  stopifnot(!is.null(harmonised_dat))

  harmonised_dat$pval.x <- pnorm(abs(harmonised_dat$beta.x)/harmonised_dat$se.x, lower.tail=FALSE)
  harmonised_dat$pval.y <- pnorm(abs(harmonised_dat$beta.y)/harmonised_dat$se.y, lower.tail=FALSE)

  # sig <- harmonised_dat[harmonised_dat$pval.x < 0.05/nrow(harmonised_dat), ]
  sig <- harmonised_dat
  sig$wr <- sig$beta.y/sig$beta.x
  sig$wr.se <- sig$se.y/abs(sig$beta.x)
  sig$biv <- NA
  sig$biv.se <- NA

  for(pop in unique(sig$pops)) {
    Qjpval <- x1$mrres$coefficients$Qjpval[x1$mrres$coefficients$pops == pop]
    if(Qjpval < 0.05) {
      b <- x1$mrres$coefficients$Estimate[x1$mrres$coefficients$pops == pop]
      se <- x1$mrres$coefficients$`Std. Error`[x1$mrres$coefficients$pops == pop]
    } else {
      b <- x1$mrres$coefficients$Estimate[x1$mrres$coefficients$pops == "All"]
      se <- x1$mrres$coefficients$`Std. Error`[x1$mrres$coefficients$pops == "All"]
    }
    sig$biv[sig$pops == pop] <- b
    sig$biv.se[sig$pops == pop] <- se
  }

  sig$dif <- b - sig$wr
  sig$dif.se <- sqrt(sig$biv.se^2 + sig$wr.se^2)
  sig <- sig %>% mutate(Qj = 1/wr.se^2 * (biv - wr)^2, Qjpval = pchisq(Qj, 1, lower.tail=FALSE))
  keepsnp <- sig[p.adjust(sig$Qjpval, "fdr") < 0.2,]$SNP %>% unique
  

  sig2 <- sig %>% dplyr::filter(SNP %in% keepsnp)
  het_of_outliers <- dplyr::group_by(sig2, SNP) %>% 
    dplyr::do({
      a <- fixed_effects_meta_analysis(.$dif, .$dif.se, 8)
      tibble(
        SNP=.$SNP,
        pop=.$pops,
        dif=.$dif,
        dif.se=.$wr.se,
        meandif=a$beta,
        meandif.se=a$se,
        Q=a$Q,
        Qpval=a$Qpval,
        Qj=a$Qj,
        Qjpval=a$Qjpval
      )
    })

  # sig2 %>% dplyr::select(SNP, pops, dif) %>% tidyr::pivot_wider(names_from=pops, values_from=dif) %>% dplyr::select(!SNP) %>% pairs

  ove <- lapply(unique(sig2$pops), \(pop) {
    s <- subset(sig2, pops == pop & Qjpval < 0.05)
    if(nrow(s) == 0) return(NULL)
    lapply(unique(sig2$pops), \(pop2) {
      if(pop == pop2) return(NULL)
      sig3 <- subset(sig2, SNP %in% s$SNP & pops==pop2) %>%
        dplyr::inner_join(s, ., by="SNP") %>% {
          prop_overlap(b_disc=.$dif.x, b_rep=.$dif.y, se_disc=.$dif.se.x, se_rep=.$dif.se.y, 0.05)
        }
      return(sig3$res %>% mutate(disc=pop, rep=pop2))
    }) %>% bind_rows()
  }) %>% bind_rows() %>% dplyr::select(disc, rep, dplyr::everything())

  self$pleiotropy_agreement <- ove
  self$pleiotropy_outliers <- sig2
  self$pleiotropy_Q_outliers <- het_of_outliers

})


CAMERA$set("public", "plot_pleiotropy", function(dat = self$pleiotropy_outliers) {
  combs <- expand.grid(pop1=unique(dat$pops), pop2=unique(dat$pops), stringsAsFactors=FALSE) %>% dplyr::filter(pop1 > pop2)
  combs <- dplyr::inner_join(combs, dat %>% dplyr::select(SNP, pop1=pops, dif1=dif, dif1.se=dif.se), relationship = "many-to-many") %>% dplyr::inner_join(., dat %>% dplyr::select(SNP, pop2=pops, dif2=dif, dif2.se=dif.se), relationship = "many-to-many")

  ggplot2::ggplot(combs %>% dplyr::filter(abs(dif1) < 10 & abs(dif2) < 10), ggplot2::aes(x=dif1, y=dif2)) +
  ggplot2::geom_point() +
  ggplot2::geom_vline(xintercept=0, linetype="dotted") +
  ggplot2::geom_hline(yintercept=0, linetype="dotted") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=dif2-1.96*dif2.se, ymax=dif2+1.96*dif2.se), width=0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin=dif1-1.96*dif1.se, xmax=dif1+1.96*dif1.se), height=0) +
  facet_grid(pop2 ~ pop1) +
  ggplot2::geom_smooth(method="lm") +
  ggplot2::labs(x="Deviation from MR estimate (pop 1)", y="Deviation from MR estimate (pop 2)")

  # combs %>% group_by(pop1, pop2) %>%
  #   dplyr::do({
  #     a <- summary(lm(dif1 ~ dif2, ., weight=1/(dif1.se^2 + dif2.se^2)))$coef
  #     tibble(b=a[2,1], se=a[2,2], pval=a[2,4])
  #   })
}) 

CAMERA$set("public", "plot_pleiotropy_heterogeneity", function(dat = self$pleiotropy_Q_outliers, pthresh=0.05) {
  dat <- subset(dat, Qpval < pthresh)
  if(nrow(dat) == 0) {
    message("No SNPs with outlier Q pval < ", pthresh)
    return(NULL)
  }
  dat %>%
    ggplot2::ggplot(., ggplot2::aes(x=dif, y=pop)) +
    ggplot2::geom_vline(xintercept=0, linetype="dotted") +
    ggplot2::geom_point(ggplot2::aes(colour=pop)) +
    ggplot2::geom_errorbarh(ggplot2::aes(colour=pop, xmin=dif-1.96*dif.se, xmax=dif+1.96*dif.se), height=0) +
    ggplot2::facet_wrap(~ SNP, scale="free_x") +
    ggplot2::labs(x="Deviation from MR estimate", y="Population")
}) 















######################















#' @description
#' Return a list of outlying SNPs in each population
#' @param harmonised_dat Harmonised dataset obtained by using \code{harmonised_dat()}
#' @return List of the pleiotropic SNPs available in x$pleiotropic_snps. Plot for the distribution of the pleiotropic SNPs based on a data frame in x$pleiotropy_dat.
CAMERA$set("public", "pleiotropy_deprecated", function(harmonised_dat = self$harmonised_dat_sem) {
  stopifnot(!is.null(harmonised_dat))

  # harmonised_dat$pval.x <- pnorm(abs(harmonised_dat$beta.x)/harmonised_dat$se.x, lower.tail=FALSE)
  # harmonised_dat$pval.y <- pnorm(abs(harmonised_dat$beta.y)/harmonised_dat$se.y, lower.tail=FALSE)

  sig <- harmonised_dat[harmonised_dat$p1 < 5 * 10^-8, ]

  d1 <- RadialMR::format_radial(BXG = sig$x1, BYG = sig$y1, seBXG = sig$xse1, seBYG = sig$yse1, RSID = sig$SNP)
  invisible(capture.output(o1 <- RadialMR::ivw_radial(d1, alpha = 1, weights = 3)))
  o1$outliers$residuals <- summary(lm(BetaWj ~ -1 + Wj, o1$data))$residuals

  d2 <- RadialMR::format_radial(BXG = sig$x2, BYG = sig$y2, seBXG = sig$xse2, seBYG = sig$yse2, RSID = sig$SNP)
  invisible(capture.output(o2 <- RadialMR::ivw_radial(d2, alpha = 1, weights = 3)))
  o2$outliers$residuals <- summary(lm(BetaWj ~ -1 + Wj, o2$data))$residuals

  dat <- merge(o1$outliers, o2$outliers, by = "SNP") %>%
    dplyr::mutate(
      sigx = p.adjust(p.value.x, "fdr") < 0.5,
      sigy = p.adjust(p.value.y, "fdr") < 0.5,
      outlier = dplyr::case_when(
        sigx & sigy ~ "Both",
        sigx & !sigy ~ "pop1",
        !sigx & sigy ~ "pop2",
        !sigx & !sigy ~ "None"
      )
    ) %>%
    dplyr::mutate(outp1 = dplyr::if_else(outlier == "pop1", 1, dplyr::if_else(outlier == "both", 1, 0))) %>%
    dplyr::mutate(outp2 = dplyr::if_else(outlier == "pop2", 1, dplyr::if_else(outlier == "both", 1, 0)))

  invisible(self$pleiotropy_dat <- dat %>% dplyr::select("SNP", "outlier", "outp1", "outp2"))

  out <- list()
  out$both <- subset(dat$SNP, dat$outlier == "Both")
  out$pop1 <- subset(dat$SNP, dat$outlier == "pop1")
  out$pop2 <- subset(dat$SNP, dat$outlier == "pop2")
  out$none <- subset(dat$SNP, dat$outlier == "None")

  invisible(self$pleiotropic_snps <- out)

  dat %>%
    ggplot2::ggplot(., ggplot2::aes(x = Q_statistic.x, y = Q_statistic.y)) +
    ggplot2::geom_point(ggplot2::aes(colour = outlier)) +
    # ggplot2::geom_smooth(method = "lm", se=FALSE, color="gray", alpha = .2, size = 0.2, formula = y ~ x) +
    ggplot2::scale_colour_brewer(type = "qual") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::ylab("Q statistics in pop2") +
    ggplot2::xlab("Q statistics in pop1")
})

#' @description
#' Estimate population specificity of pleiotric SNPs
#' @param harmonised_dat Harmonised dataset obtained by using \code{harmonised_dat()}
#' @param sem_result MR-SEM result obtained by \code{x$perform_basic_sem()}
#' @param pleioropy List of pleiotropic SNPs obtained using \code{x$pleiotropy()}
#' @return Table of the results. Summary of the results available in x$instrument_pleiotropy_summary.
CAMERA$set("public", "pleiotropy_specificity_deprecated", function(harmonised_dat = self$harmonised_dat_sem, sem_result = self$sem_result, pleioropy = self$pleiotropy_dat) {
  stopifnot(!is.null(harmonised_dat))
  stopifnot(!is.null(sem_result))

  sig <- harmonised_dat[harmonised_dat$p1 < 5 * 10^-8, ]

  if (sem_result$aic[5] - sem_result$aic[6] <= -2) {
    d <- sig %>%
      dplyr::mutate(
        wald1 = y1 / x1, wald.se1 = yse1 / abs(x1), wald2 = y2 / x2, wald.se2 = yse2 / abs(x2),
        ivw1 = sem_result$bivhat[5], ivw.se1 = sem_result$se[5], ivw2 = sem_result$bivhat[5], ivw.se2 = sem_result$se[5]
      )
  }

  if (is.na(sem_result$aic[6])) {
    message("Caution: SE is not properly estimated for model 2. The estimates from model 1 are used.")
    d <- sig %>%
      dplyr::mutate(
        wald1 = y1 / x1, wald.se1 = yse1 / abs(x1), wald2 = y2 / x2, wald.se2 = yse2 / abs(x2),
        ivw1 = sem_result$bivhat[5], ivw.se1 = sem_result$se[5], ivw2 = sem_result$bivhat[5], ivw.se2 = sem_result$se[5]
      )
  }

  if (sem_result$aic[5] - sem_result$aic[6] > -2) {
    d <- sig %>%
      dplyr::mutate(
        wald1 = y1 / x1, wald.se1 = yse1 / abs(x1), wald2 = y2 / x2, wald.se2 = yse2 / abs(x2),
        ivw1 = sem_result$bivhat[6], ivw.se1 = sem_result$se[6], ivw2 = sem_result$bivhat[7], ivw.se2 = sem_result$se[7]
      )
  }

  pop1 <- list()
  pop2 <- list()
  for (i in 1:nrow(d))
  {
    pop1[[i]] <- private$bootstrap(d$wald1[i], d$wald.se1[i], d$ivw1[i], d$ivw.se1[i], nboot = 1000)
    pop2[[i]] <- private$bootstrap(d$wald2[i], d$wald.se2[i], d$ivw2[i], d$ivw.se2[i], nboot = 1000)
  }
  pop1 <- do.call("rbind", pop1) %>%
    as.data.frame() %>%
    dplyr::select(pleio.p1 = pleio, sd.p1 = sd)
  pop2 <- do.call("rbind", pop2) %>%
    as.data.frame() %>%
    dplyr::select(pleio.p2 = pleio, sd.p2 = sd)

  d <- cbind(d, pop1, pop2) %>%
    merge(., pleioropy, by = "SNP")

  p1 <- d %>%
    dplyr::select("SNP", x = "x1", xse = "xse1", p = "p1", y = "y1", yse = "yse1", wald = "wald1", wald.se1 = "wald.se1", pleio = "pleio.p1", sd = "sd.p1", out = "outp1") %>%
    dplyr::mutate(id = self$exposure_ids[1])
  p2 <- d %>%
    dplyr::select("SNP", x = "x2", xse = "xse2", p = "p2", y = "y2", yse = "yse2", wald = "wald2", wald.se1 = "wald.se2", pleio = "pleio.p2", sd = "sd.p2", out = "outp2") %>%
    dplyr::mutate(id = self$exposure_ids[2])

  mer <- rbind(p1, p2)

  o <- lapply(self$exposure_ids, function(i) {
    m <- subset(mer, id == i & out == 1)
    other_ids <- self$exposure_ids[!self$exposure_ids %in% i]
    o <- lapply(other_ids, function(j) {
      n <- subset(mer, id == j & SNP %in% m$SNP) %>% dplyr::arrange(SNP)
      o <- n %>%
        {
          prop_overlap(m$pleio, .$pleio, sqrt(m$sd), sqrt(.$sd), alpha = 0.05)
        }
      o$res <- o$res %>%
        dplyr::mutate(discovery = i, replication = j) %>%
        dplyr::select(discovery, replication, dplyr::everything())
      o$variants <- o$variants %>%
        dplyr::mutate(
          discovery = i,
          replication = j,
          rsid = n$rsid,
          p = n$p,
          distinct = sig > 0.8 & p > 0.1
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

  res <- lapply(o, function(x) x$overall) %>% dplyr::bind_rows()

  pval <- list()
  for (i in 1:(nrow(res) / 2))
  {
    if (res$nsnp[c(FALSE, TRUE)][i] != 0) {
      pval[i] <- round(binom.test(res$value[c(FALSE, TRUE)][i], res$nsnp[c(FALSE, TRUE)][i], p = res$value[c(TRUE, FALSE)][i] / res$nsnp[c(FALSE, TRUE)][i])$p.value, 3)
    }
    if (res$nsnp[c(FALSE, TRUE)][i] == 0) {
      pval[i] <- NA
    }
  }
  pval <- tibble::tibble(pvalue_diff = unlist(pval, use.names = FALSE))
  pval <- pval[rep(1:nrow(pval), each = 2), ]

  self$instrument_pleiotropy_summary <- cbind(res, pval) %>% dplyr::mutate(pvalue_diff = ifelse(dplyr::row_number() %% 2, pvalue_diff, "."))

  print(self$instrument_pleiotropy_summary)
})
