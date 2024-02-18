#' @description
#' The function standardise the betas and SEs for the instruments-exposure/outcome associations when unit informaiton is not matched across the populations. In case of large differences in genetic effects observed between the populations (e.g. due to sample size difference), the function scales the betas and SEs.
#' @param dat Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param standardise_unit Use this option if unit information is not matched.
#' @param standardise_scale Use this option if genetic effects are substantially different due to study power.
#' @param scaling_method Choose the methods to obtain scaling units (MR estimates of exposure 1 and exposure 2 or outcome 1 and outcome 2). Default is "simple_mode".
#' @return Data frame in x$standardised_instrument_raw; x$standardised_instrument_maxz; x$standardised_instrument_susie; x$standardised_instrument_paintor; x$standardised_instrument_mscaviar; x$standardised_outcome
CAMERA$set("public", "standardise_data", function(dat = self$instrument_raw, standardise_unit = FALSE, standardise_scale = FALSE, scaling_method = "simple_mode") {
  if (standardise_unit == TRUE) {
    if (!any(names(dat) %in% c("beta.outcome"))) {
      exp <- dat
      d <- exp %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(units = units[1])

      if (any(is.na(exp$eaf))) {
        exp <- private$allele_frequency(dat = exp)
      }

      if (!any(d$units %in% c("log odds"))) {
        exp <- exp %>%
          dplyr::group_by(id) %>%
          dplyr::mutate(units = dplyr::na_if(units, "NA")) %>%
          dplyr::mutate(units = replace(units, is.na(units), "temp")) %>%
          dplyr::mutate(estimated_sd = mean(TwoSampleMR::estimate_trait_sd(beta, se, samplesize, eaf), na.rm = TRUE)) %>%
          dplyr::mutate(estimated_sd = replace(estimated_sd, units == "SD", 1))
      }

      if (!any(is.na(exp$estimated_sd))) {
        stopifnot(!any(is.na(exp$estimated_sd)))
        exp$beta <- exp$beta / exp$estimated_sd
        exp$se <- exp$se / exp$estimated_sd
        exp$units <- "SD"

        exp <- exp %>% as.data.frame()

        if (any(exp$method[[1]] %in% c("raw"))) {
          self$standardised_instrument_raw <- exp
        }
        if (any(exp$method[[1]] %in% c("maxz"))) {
          self$standardised_instrument_maxz <- exp
        }
        if (any(exp$method[[1]] %in% c("susie"))) {
          self$standardised_instrument_susie <- exp
        }
        if (any(exp$method[[1]] %in% c("paintor"))) {
          self$standardised_instrument_paintor <- exp
        }
        if (any(exp$method[[1]] %in% c("mscaviar"))) {
          self$standardised_instrument_mscaviar <- exp
        }
      }
    }

    if (any(names(dat) %in% c("beta.outcome"))) {
      out <- dat
      d <- out %>%
        dplyr::group_by(id.outcome) %>%
        dplyr::summarise(units = units.outcome[1])

      if (any(is.na(out$eaf.outcome))) {
        out <- private$allele_frequency(dat = out)
      }

      if (!any(d$units %in% c("log odds"))) {
        out <- out %>%
          dplyr::group_by(id.outcome) %>%
          dplyr::mutate(units.outcome = dplyr::na_if(units.outcome, "NA")) %>%
          dplyr::mutate(units.outcome = replace(units.outcome, is.na(units.outcome), "temp")) %>%
          dplyr::mutate(estimated_sd = mean(TwoSampleMR::estimate_trait_sd(beta.outcome, se.outcome, samplesize.outcome, eaf.outcome), na.rm = TRUE)) %>%
          dplyr::mutate(estimated_sd = replace(estimated_sd, units.outcome == "SD", 1))
      }

      if (any(is.na(out$estimated_sd))) {
        stopifnot(!any(is.na(out$estimated_sd)))
        out$beta.outcome <- out$beta.outcome / out$estimated_sd
        out$se.outcome <- out$se.outcome / out$estimated_sd
        out$units.outcome <- "SD"
      }

      self$standardised_outcome <- out %>% as.data.frame()
    }
  }

  if (standardise_scale == TRUE) {
    if (!any(names(dat) %in% c("beta.outcome"))) {
      stopifnot(!is.null(self$instrument_maxz))
      if (!length(self$instrument_maxz$units) == 1) {
        invisible(capture.output(self$standardise_data(dat = self$instrument_maxz, standardise_unit = TRUE, standardise_scale = FALSE)))
        inst <- self$standardised_instrument_maxz
        invisible(capture.output(scale <- self$instrument_heterogeneity(instrument = inst, method = scaling_method, outlier_removal = TRUE)))
      }
      if (length(self$instrument_maxz$units) == 1) {
        invisible(capture.output(scale <- self$instrument_heterogeneity(instrument = self$instrument_maxz, method = scaling_method, outlier_removal = TRUE)))
      }
      bxx <- scale$agreement[1]
      if (!is.null(self[[paste0("standardised_instrument_", dat$method[[1]])]])) {
        dat <- self[[paste0("standardised_instrument_", dat$method[[1]])]]
      }
      exp <- dat %>%
        dplyr::mutate(original_beta = beta) %>%
        dplyr::mutate(original_se = se) %>%
        dplyr::mutate(beta = dplyr::case_when(id == self$exposure_ids[1] ~ beta * bxx, TRUE ~ beta)) %>%
        dplyr::mutate(se = dplyr::case_when(id == self$exposure_ids[1] ~ se * bxx, TRUE ~ se)) %>%
        as.data.frame()

      if (any(exp$method[[1]] %in% c("raw"))) {
        self$standardised_instrument_raw <- exp
      }
      if (any(exp$method[[1]] %in% c("maxz"))) {
        self$standardised_instrument_maxz <- exp
      }
      if (any(exp$method[[1]] %in% c("susie"))) {
        self$standardised_instrument_susie <- exp
      }
      if (any(exp$method[[1]] %in% c("paintor"))) {
        self$standardised_instrument_paintor <- exp
      }
      if (any(exp$method[[1]] %in% c("mscaviar"))) {
        self$standardised_instrument_mscaviar <- exp
      }
    }

    if (any(names(dat) %in% c("beta.outcome"))) {
      stopifnot(!is.null(self$instrument_outcome))

      if (!is.null(self$standardised_outcome)) {
        dat <- self$standardised_outcome
      }

      oexp <- self$exposure_ids
      oraw <- self$instrument_raw
      omaxz <- self$instrument_maxz
      ore <- self$instrument_regions
      oz <- self$instrument_region_zscores

      self$exposure_ids <- self$outcome_ids
      suppressMessages(self$extract_instruments(exposure_ids = self$exposure_ids))
      suppressMessages(self$extract_instrument_regions(instrument_raw = self$instrument_raw, exposure_ids = self$exposure_ids))
      suppressMessages(self$scan_regional_instruments(instrument_raw = self$instrument_raw, instrument_regions = self$instrument_regions))
      invisible(capture.output(scale <- self$instrument_heterogeneity(instrument = self$instrument_maxz, method = scaling_method)))

      if (is.na(scale$agreement[1])) {
        stop("Scaling units can't be estimated")
      }

      byy <- scale$agreement[1]
      out <- dat
      out <- out %>%
        dplyr::mutate(original_beta = beta.outcome) %>%
        dplyr::mutate(original_se = se.outcome) %>%
        dplyr::mutate(beta.outcome = dplyr::case_when(id.outcome == self$outcome_ids[1] ~ beta.outcome * byy, TRUE ~ beta.outcome)) %>%
        dplyr::mutate(se.outcome = dplyr::case_when(id.outcome == self$outcome_ids[1] ~ se.outcome * byy, TRUE ~ se.outcome)) %>%
        as.data.frame()

      self$exposure_ids <- oexp
      self$instrument_raw <- oraw
      self$instrument_maxz <- omaxz
      self$instrument_regions <- ore
      self$instrument_region_zscores <- oz

      self$standardised_outcome <- out
    }
  }
  invisible(self)
})

#' @importFrom dplyr group_by mutate
#' @importFrom tidyr replace_na
CAMERA$set("private", "sd_standardise", function(dat = dat) {
  d <- dat %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(units = units[1])
  if (!any(d$units %in% c("log odds"))) {
    dat <- dat %>%
      dplyr::group_by(id) %>%
      tidyr::replace_na(list(units = "temp")) %>%
      dplyr::mutate(estimated_sd = mean(TwoSampleMR::estimate_trait_sd(beta, se, samplesize, eaf), na.rm = TRUE)) %>%
      dplyr::mutate(estimated_sd = replace(estimated_sd, units == "SD", 1))
  }
  if (any(!is.na(dat$estimated_sd))) {
    stopifnot(!any(is.na(dat$estimated_sd)))
    dat$beta <- dat$beta / dat$estimated_sd
    dat$se <- dat$se / dat$estimated_sd
    dat$units <- "SD"
  }
})


#' @importFrom dplyr mutate group_by bind_rows
#' @importFrom ieugwasr afl2_rsid
CAMERA$set("private", "allele_frequency", function(dat = dat) {
  if (!any(names(dat) %in% c("beta.outcome"))) {
    dat <- dat %>%
      dplyr::mutate(pops = factor(ifelse(id == self$exposure_ids[[1]], self$pops[[1]], self$pops[[2]])))
    index <- which(is.na(dat$eaf2))
    id <- unique(dat$id[index])
    pop <- unique(dat$pops[index])

    if (length(id) == 1) {
      dat <- dat %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(eaf = ifelse(is.na(eaf), ieugwasr::afl2_rsid(rsid)[[paste0("AF.", pop)]], eaf)) %>%
        as.data.frame()
    }

    if (length(id) > 1) {
      af <- list()
      for (i in pop) {
        af[[i]] <- dat %>%
          subset(., pops == i) %>%
          dplyr::mutate(eaf = ifelse(is.na(eaf), ieugwasr::afl2_rsid(rsid)[[paste0("AF.", pop)]], eaf)) %>%
          as.data.frame()
      }
      dat <- af %>% dplyr::bind_rows()
    }
  }

  if (any(names(dat) %in% c("beta.outcome"))) {
    dat <- dat %>%
      dplyr::mutate(pops = factor(ifelse(id.outcome == self$outcome_ids[[1]], self$pops[[1]], self$pops[[2]])))

    index <- which(is.na(dat$eaf.outcome))
    id <- unique(dat$id.outcome[index])
    pop <- unique(dat$pops[index])

    if (length(id) == 1) {
      dat <- dat %>%
        dplyr::mutate(eaf.outcome = ifelse(is.na(eaf.outcome), ieugwasr::afl2_rsid(SNP)[[paste0("AF.", pop)]], eaf.outcome)) %>%
        as.data.frame()
    }

    if (length(id) > 1) {
      af <- list()
      for (i in pop) {
        af[[i]] <- dat %>%
          subset(., pops == i) %>%
          dplyr::mutate(eaf.outcome = ifelse(is.na(eaf.outcome), ieugwasr::afl2_rsid(SNP)[[paste0("AF.", i)]], eaf.outcome)) %>%
          as.data.frame()
      }
      dat <- af %>% dplyr::bind_rows()
    }
  }
  return(dat)
})
