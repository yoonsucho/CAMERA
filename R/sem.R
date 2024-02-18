

#' @description
#' Perform basic SEM analysis of the joint estimates in multiple ancestries.
#' @param harmonised_dat Harmonised dataset obtained by using \code{harmonised_dat()}
#' @return Table of the SEM results. Summary of the results available in x$sem_result.
CAMERA$set("public", "perform_basic_sem", function(harmonised_dat = self$harmonised_dat_sem) {
  d <- harmonised_dat %>%
    dplyr::mutate(r1 = y1 / x1) %>%
    dplyr::mutate(r2 = y2 / x2) %>%
    dplyr::mutate(w1 = sqrt(x1^2 / yse1^2)) %>%
    dplyr::mutate(w2 = sqrt(x2^2 / yse2^2)) %>%
    dplyr::mutate(o1 = r1 * w1) %>%
    dplyr::mutate(o2 = r2 * w2)

  out <- list()
  out$ivw1 <- TwoSampleMR::mr_ivw(d$x1, d$y1, d$xse1, d$yse1) %>%
    {tibble::tibble(Methods = "IVW", pop = "1", nsnp = nrow(d), bivhat = .$b, se = .$se, pval = .$pval)}
  out$ivw2 <- TwoSampleMR::mr_ivw(d$x2, d$y2, d$xse2, d$yse2) %>%
    {tibble::tibble(Methods = "IVW", pop = "2", nsnp = nrow(d), bivhat = .$b, se = .$se, pval = .$pval)}
  out$rm1 <- summary(lm(o1 ~ -1 + w1, data = d)) %>%
    {tibble::tibble(Methods = "RadialIVW", pop = "1", nsnp = nrow(d), bivhat = .$coef[1, 1], se = .$coef[1, 2], pval = .$coef[1, 4])}
  out$rm2 <- summary(lm(o2 ~ -1 + w2, data = d)) %>%
    {tibble::tibble(Methods = "RadialIVW", pop = "2", nsnp = nrow(d), bivhat = .$coef[1, 1], se = .$coef[1, 2], pval = .$coef[1, 4])}

  out$semA <- self$runsem("
                                y1 ~ biv*x1
                                y2 ~ biv*x2
                                ", d, "UnweightedSEMa")[1, ] %>%
    dplyr::mutate(pop = replace(pop, pop == 1, "1=2"))

  out$semB <- self$runsem("
                                y1 ~ biv_1*x1
                                y2 ~ biv_2*x2
                                ", d, "UnweightedSEMb")

  out$modC <- self$runsem("
                                o1 ~ biv*w1
                                o2 ~ biv*w2
                                ", d, "RadialSEMa")[1, ] %>%
    dplyr::mutate(pop = replace(pop, pop == 1, "1=2"))

  out$modD <- self$runsem("
                                o1 ~ biv_1*w1
                                o2 ~ biv_2*w2
                                ", d, "RadialSEMb")

  invisible(self$sem_result <- out %>% dplyr::bind_rows())
  print(self$sem_result)
})


#' @importFrom tibble tibble
#' @importFrom dplyr mutate
CAMERA$set("public", "runsem", function(model, data, modname) {
  mod <- lavaan::sem(model, data = data)
  invisible(capture.output(mod <- lavaan::summary(mod, fit.measures = TRUE)))
  o <- tibble::tibble(
    Methods = modname,
    pop = 1:2,
    nsnp = nrow(data),
    bivhat = mod$pe$est[1:2],
    se = mod$pe$se[1:2],
    pval = mod$pe$pvalue[1:2],
    aic = mod$fit["aic"]
  ) %>% dplyr::mutate(pop = as.character(pop))
  temp <<- mod
  print(mod)
  print(o)
  if (any(is.na(o$se))) {
    message("WARNING: The model convergence was not successful. No constraints were used.")
    mod <- lavaan::sem(model, data = data, check.gradient = FALSE)
    invisible(capture.output(mod <- lavaan::summary(mod, fit.measures = TRUE)))
    o <- tibble::tibble(
      Methods = modname,
      pop = 1:2,
      nsnp = nrow(data),
      bivhat = mod$pe$est[1:2],
      se = mod$pe$se[1:2],
      pval = mod$pe$pvalue[1:2],
      aic = mod$fit["aic"]
    ) %>% dplyr::mutate(pop = as.character(pop))
  }
  return(o)
})


CAMERA$set("private", "jackknife2", function(x, theta, ...) {
  call <- match.call()
  n <- length(x)
  u <- rep(0, n)
  for (i in 1:n) {
    u[i] <- theta(x[-i], ...)
  }
  thetahat <- theta(x, ...)
  jack.bias <- (n - 1) * (mean(u) - thetahat)
  jack.se <- sqrt(((n - 1) / n) * sum((u - mean(u))^2))
  jack.sd <- sqrt((1 / n) * sum((u - mean(u))^2))
  return(list(
    jack.se = jack.se, jack.sd = jack.sd, jack.bias = jack.bias, jack.values = u,
    call = call
  ))
})


#' @importFrom tibble tibble
CAMERA$set("private", "bootstrap_diff", function(nboot, slope, slope_se, b_out, b_out_se, b_exp, b_exp_se) {
  expected_b_out <- b_exp * slope
  diff <- b_out - expected_b_out

  # bootstrap to get diff_se
  boots <- tibble::tibble(
    B_EXP = rnorm(nboot, mean = b_exp, sd = b_exp_se),
    B_OUT = rnorm(nboot, mean = b_out, sd = b_out_se),
    SLOPE = rnorm(nboot, mean = slope, sd = slope_se),
    DIFF = B_OUT - B_EXP * SLOPE
  )
  diff_se <- sd(boots$DIFF)
  return(list(
    diff = diff,
    diff_se = diff_se
  ))
})


CAMERA$set("private", "bootstrap", function(wr, wr.se, ivw, ivw.se, nboot = 1000) {
  res <- rnorm(nboot, wr, wr.se) - rnorm(nboot, ivw, ivw.se)
  pe <- wr - ivw
  return(c(pleio = pe, sd = sd(res)))
})
