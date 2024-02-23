# Once a set of instruments is defined then extract the outcome data
# Could use
# - raw results from extract_instruments
# - maximised associations from extract_instrument_regions
# - finemapped hits from susie_finemap_regions
# - finemapped hits from mscaviar_finemap_regions
#' @description
#' The function extracts summary statistics of given a list of instruments and the outcomes
#' @param exp Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param p_exp Statistical threshold to determine significance. Default is "bonferroni", which is eqaul to 0.05/number of the instruments.
#' @return Data frame in x$instrument_outcome
CAMERA$set("public", "make_outcome_data", function(exp = self$instrument_raw, p_exp = 0.05 / nrow(exp)) {
  rsids <- unique(exp$rsido)
  if(!is.null(self$instrument_outcome)) {
    message("Appending new outcome data to existing outcome data")
    rsids <- rsids[!rsids %in% self$instrument_outcome$SNP]
    message(length(rsids))
  }
  out <- TwoSampleMR::extract_outcome_data(snps = rsids, outcomes = self$outcome_ids) %>%
    generate_vid(., ea = "effect_allele.outcome", nea = "other_allele.outcome", eaf = "eaf.outcome", beta = "beta.outcome", rsid = "SNP", chr = "chr", position = "pos")
  print(str(out))
  suppressMessages(out <- TwoSampleMR::add_metadata(out, cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")) %>%
    # dplyr::select(rsid=SNP, chr, position=pos, id=id.outcome, beta=beta.outcome, se=se.outcome, p=pval.exposure, ea=effect_allele.outcome, nea=other_allele.outcome, eaf=eaf.outcome, units=units.outcome, samplesize=contains("size")) %>%
    as.data.frame())
  if(!is.null(self$instrument_outcome)) {
    out_orig <- self$instrument_outcome
    out <- dplyr::bind_rows(out_orig, out) %>% subset(!duplicated(paste(SNP, id.outcome)))
  }
  self$instrument_outcome <- out %>% dplyr::arrange(., chr, SNP)
  invisible(self)
})


# Extract a set of instruments from locally-supplied outcome regional data
# Could use
# - raw results from extract_instruments
# - maximised associations from extract_instrument_regions
# - finemapped hits from susie_finemap_regions
# - finemapped hits from mscaviar_finemap_regions
#' @description
#' The function extracts summary statistics of given a list of instruments and the outcomes
#' @param exp Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_fema, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param p_exp Statistical threshold to determine significance. Default is "bonferroni", which is eqaul to 0.05/number of the instruments.
#' @return Data frame in x$instrument_outcome
CAMERA$set("public", "make_outcome_local", function(exp = self$instrument_raw, out = self$instrument_outcome_regions, p_exp = 0.05 / nreow(exp)) {
  out <- lapply(out, \(x) {
    lapply(x, \(y) {
      subset(y, rsid %in% exp$rsid)
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  self$instrument_outcome <- out %>% dplyr::arrange(., chr, rsid)
})


#' Harmonise exposure and outcome datasets
#' @description
#' This function harmonises the alleles and effects between the exposure and outcome.
#' @param exp Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_fema, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param out Intsruments for the outcome by using \code{make_outcome_data()}. Default is x$instrument_outcome.
#' @return Data frame in x$harmonised_dat_sem
CAMERA$set("public", "harmonise", function(exp = self$instrument_raw, out = self$instrument_outcome) {
  self$set_summary()
  if(self$source == "OpenGWAS") {
    exp <- dplyr::left_join(exp, subset(self$summary, select=c(exposure_ids, pops)), by=c("id"="exposure_ids")) %>%
      dplyr::select(SNP = rsid, pops, beta, se)
    out <- dplyr::left_join(out, subset(self$summary, select=c(outcome_ids, pops)), by=c("id.outcome" = "outcome_ids")) %>%
      dplyr::select(SNP, pops, beta=beta.outcome, se=se.outcome)
    print(str(exp))
    print(str(out))
    dat <- dplyr::inner_join(exp, out, by=c("pops", "SNP"))
  } else {
    dat <- dplyr::inner_join(exp, out, by=c("pop", "rsid")) %>% dplyr::rename(pops="pop", SNP="rsid")
  }
  self$harmonised_dat <- dat
  print(str(dat))
})

#' Generate summary of exposure, outcome and population metadata
#' 
#' @return data frame
CAMERA$set("public", "set_summary", function() {
  self$summary <- dplyr::tibble(
    pops = self$pops,
    exposure_ids = self$exposure_ids,
    outcome_ids = self$outcome_ids,
    source=self$source
  )
  print(self$summary)
})
