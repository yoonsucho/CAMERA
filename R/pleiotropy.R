#' Estimate similarity of horizontal pleiotropy across ancestries
#' 
#' For each ancestry, identify outliers in the MR analysis based on per-variatn Q statistics. Then estimate the deviation from the main estimates for all outliers across all ancestries. Finally, determine if the pleiotropy deviation is consistent across all ancestries
#' 
#' @param harmonised_dat Outcome from `harmonise` function
#' @param mrres Outcome from `cross_estimate`
#' 
#' @return 
#'    - data frame of outliers `pleiotropy_outliers`
#'    - data frame of heterogeneity for each outlier / ancestry combination `pleiotropy_Q_outliers`
#'    - data frame of agreement of outlier effects `pleiotropy_agreement`
CAMERA$set("public", "pleiotropy", function(harmonised_dat = self$harmonised_dat, mrres=self$mrres) {
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
    Qjpval <- mrres$Qjpval[mrres$pops == pop]
    if(Qjpval < 0.05) {
      b <- mrres$Estimate[mrres$pops == pop]
      se <- mrres$`Std. Error`[mrres$pops == pop]
    } else {
      b <- mrres$Estimate[mrres$pops == "All"]
      se <- mrres$`Std. Error`[mrres$pops == "All"]
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
        dplyr::inner_join(s, ., by="SNP") 
        if(nrow(sig3) == 0) return(NULL)
        sig3 <- sig3 %>% {
          prop_overlap(b_disc=.$dif.x, b_rep=.$dif.y, se_disc=.$dif.se.x, se_rep=.$dif.se.y, 0.05)
        }
      return(sig3$res %>% mutate(disc=pop, rep=pop2))
    }) %>% bind_rows()
  }) %>% bind_rows() %>% dplyr::select(disc, rep, dplyr::everything())

  self$pleiotropy_agreement <- ove
  self$pleiotropy_outliers <- sig2
  self$pleiotropy_Q_outliers <- het_of_outliers

})

#' Plot pleiotropy results
#' 
#' @param dat Output from `pleiotropy` - `pleiotropy_outliers`
#' 
#' @return plot
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

#' Plot pleiotropy results per variant
#' 
#' @param dat Output from `pleiotropy` - `pleiotropy_Q_outliers`
#' @param pthresh p-value from Q statistic for inclusion in plots
#' 
#' @return plot
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

