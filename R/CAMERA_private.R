#' @importFrom tibble tibble
CAMERA$set("private", "prop_overlap", function(b_disc, b_rep, se_disc, se_rep, alpha)
{
  p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))
  p_rep <- pnorm(abs(b_rep)/se_rep, lower.tail=FALSE)
  res <- tibble::tibble(
    nsnp=length(b_disc),
    metric=c("Sign", "Sign", "P-value", "P-value"),
    datum=c("Expected", "Observed", "Expected", "Observed"),
    value=c(sum(p_sign, na.rm=TRUE), sum(sign(b_disc) == sign(b_rep)), sum(p_sig, na.rm=TRUE), sum(p_rep < alpha, na.rm=TRUE))
  )
  return(list(res=res, variants=dplyr::tibble(sig=p_sig, sign=p_sign)))
})

#' @importFrom dplyr bind_rows mutate arrange
#' @importFrom tibble tibble
CAMERA$set("private", "susie_overlaps", function(su1, su2)
{
  l <- list()
  k <- 1
  s1 <- su1$sets$cs
  s2 <- su2$sets$cs
  for(i in 1:length(s1))
  {
    for(j in 1:length(s2))
    {
      if(any(s1[[i]] %in% s2[[j]]))
      {
        ind <- s1[[i]] %in% s2[[j]]
        v <- s1[[i]][ind]
        l[[k]] <- tibble::tibble(s1=i, s2=j, v = v)
        k <- k + 1
      }
    }
  }
  l <- dplyr::bind_rows(l)
  if(nrow(l) > 0)
  {
    l$pip1 <- su1$pip[l$v]
    l$pip2 <- su2$pip[l$v]
    l$piprank1 <- rank(-su1$pip)[l$v] / length(su1$pip)
    l$piprank2 <- rank(-su2$pip)[l$v] / length(su2$pip)
    l <- l %>%
      dplyr::mutate(piprank=piprank1 + piprank2) %>%
      dplyr::arrange(piprank)
  }
  return(l)
})

#' @importFrom tibble tibble
#' @importFrom lavaan sem
#' @importFrom dplyr mutate
CAMERA$set("private", "runsem", function(model, data, modname)
{
  mod <- lavaan::sem(model, data=data)
  invisible(capture.output(mod <- lavaan::summary(mod, fit.measures=TRUE)))
  o <- tibble::tibble(
                      Methods=modname,
                      pop=1:2,
                      nsnp=nrow(data),
                      bivhat=mod$PE$est[1:2],
                      se=mod$PE$se[1:2],
                      pval=mod$PE$pvalue[1:2],
                      aic=mod$FIT['aic']
                    ) %>%  dplyr::mutate(pop = as.character(pop))

  if (is.na(o$se)){
    message("WARNING: The model convergence was not successful. No constraints were used.")
    mod <- lavaan::sem(model, data=data, check.gradient = FALSE)
    invisible(capture.output(mod <- lavaan::summary(mod, fit.measures=TRUE)))
    o <- tibble::tibble(
                      Methods=modname,
                      pop=1:2,
                      nsnp=nrow(data),
                      bivhat=mod$PE$est[1:2],
                      se=mod$PE$se[1:2],
                      pval=mod$PE$pvalue[1:2],
                      aic=mod$FIT['aic']
                    ) %>%  dplyr::mutate(pop = as.character(pop))
  }
  return(o)
})


CAMERA$set("private", "greedy_remove", function(r, thresh)
{
  diag(r) <- 0
  r <- abs(r)
  flag <- 1
  rem <- c()
  nom <- colnames(r)
  while(flag == 1)
  {
    count <- apply(r, 2, function(x) sum(x >= thresh))
    if(any(count > 0))
    {
      worst <- which.max(count)[1]
      rem <- c(rem, names(worst))
      r <- r[-worst,-worst]
    } else {
      flag <- 0
    }
  }
  return(which(nom %in% rem))
})


CAMERA$set("private", "jackknife2", function (x, theta, ...)
{
  call <- match.call()
  n <- length(x)
  u <- rep(0, n)
  for (i in 1:n) {
    u[i] <- theta(x[-i], ...)
  }
  thetahat <- theta(x, ...)
  jack.bias <- (n - 1) * (mean(u) - thetahat)
  jack.se <- sqrt(((n - 1)/n) * sum((u - mean(u))^2))
  jack.sd <- sqrt((1/n) * sum((u - mean(u))^2))
  return(list(jack.se = jack.se, jack.sd = jack.sd, jack.bias = jack.bias, jack.values = u,
              call = call))
})

#' @importFrom TwoSampleMR estimate_trait_sd
#' @importFrom dplyr group_by mutate
#' @importFrom tidyr replace_na
CAMERA$set("private", "sd_standardise", function (dat=dat)
{
 d <- dat %>%
        dplyr::group_by(id) %>% dplyr::summarise(units = units[1])
 if(!any(d$units %in% c("log odds"))){
     dat <- dat %>%
                 dplyr::group_by(id) %>%
                 tidyr::replace_na(list(units = "temp")) %>%
                 dplyr::mutate(estimated_sd = mean(TwoSampleMR::estimate_trait_sd(beta, se, samplesize, eaf), na.rm=TRUE)) %>%
                 dplyr::mutate(estimated_sd = replace(estimated_sd, units=="SD", 1))
 }
 if(any(!is.na(dat$estimated_sd)))
 {
   stopifnot(!any(is.na(dat$estimated_sd)))
   dat$beta <- dat$beta / dat$estimated_sd
   dat$se <- dat$se / dat$estimated_sd
   dat$units <- "SD"
 }
})

#' @importFrom dplyr mutate group_by bind_rows
#' @importFrom ieugwasr afl2_rsid
CAMERA$set("private", "allele_frequency", function (dat=dat)
{
 if(!any(names(dat) %in% c("beta.outcome")))
 {
    dat <- dat %>%
               dplyr::mutate(pops = factor(ifelse(id == self$exposure_ids[[1]], self$pops[[1]], self$pops[[2]])))
    index <- which(is.na(dat$eaf2))
    id <- unique(dat$id[index])
    pop <- unique(dat$pops[index])

    if(length(id)==1) {
      dat <- dat %>%
                dplyr::group_by(id) %>%
                dplyr::mutate(eaf = ifelse(is.na(eaf), ieugwasr::afl2_rsid(rsid)[[paste0("AF.", pop)]], eaf)) %>%
                as.data.frame()}

    if(length(id)>1) {
      af <- list()
      for (i in pop){
        af[[i]] <- dat %>%
        subset(., pops==i) %>%
                           dplyr::mutate(eaf = ifelse(is.na(eaf), ieugwasr::afl2_rsid(rsid)[[paste0("AF.", pop)]], eaf)) %>%
                           as.data.frame()}
      dat <- af %>% dplyr::bind_rows()}
  }

 if(any(names(dat) %in% c("beta.outcome")))
 {
    dat <- dat %>%
      dplyr::mutate(pops = factor(ifelse(id.outcome == self$outcome_ids[[1]], self$pops[[1]], self$pops[[2]])))

    index <- which(is.na(dat$eaf.outcome))
    id <- unique(dat$id.outcome[index])
    pop <- unique(dat$pops[index])

    if(length(id)==1) {
        dat <- dat %>%
               dplyr::mutate(eaf.outcome = ifelse(is.na(eaf.outcome), ieugwasr::afl2_rsid(SNP)[[paste0("AF.", pop)]], eaf.outcome)) %>%
               as.data.frame()}

    if(length(id)>1) {
      af <- list()
      for (i in pop){
        af[[i]] <- dat %>%
                      subset(., pops==i)  %>%
                      dplyr::mutate(eaf.outcome = ifelse(is.na(eaf.outcome), ieugwasr::afl2_rsid(SNP)[[paste0("AF.", i)]], eaf.outcome)) %>%
                      as.data.frame()}
      dat <- af %>% dplyr::bind_rows()
     }
  }
  return(dat)
})


#' @importFrom tibble tibble
CAMERA$set("private", "bootstrap_diff", function(nboot, slope, slope_se, b_out, b_out_se, b_exp, b_exp_se)
{
  expected_b_out <- b_exp * slope
  diff <- b_out - expected_b_out

  # bootstrap to get diff_se
  boots <- tibble::tibble(
    B_EXP = rnorm(nboot, mean=b_exp, sd=b_exp_se),
    B_OUT = rnorm(nboot, mean=b_out, sd=b_out_se),
    SLOPE = rnorm(nboot, mean=slope, sd=slope_se),
    DIFF = B_OUT - B_EXP * SLOPE
  )
  diff_se <- sd(boots$DIFF)
  return(list(
    diff = diff,
    diff_se = diff_se
  ))
})


CAMERA$set("private", "bootstrap", function(wr, wr.se, ivw, ivw.se, nboot=1000)
{
  res <- rnorm(nboot, wr, wr.se) - rnorm(nboot, ivw, ivw.se)
  pe <- wr - ivw
  return(c(pleio=pe, sd=sd(res)))
})
