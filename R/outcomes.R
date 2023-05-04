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
  out <- TwoSampleMR::extract_outcome_data(snps = rsids, outcomes = self$outcome_ids) %>%
    generate_vid(., ea = "effect_allele.outcome", nea = "other_allele.outcome", eaf = "eaf.outcome", beta = "beta.outcome", rsid = "SNP", chr = "chr", position = "pos")
  print(head(out))
  suppressMessages(out <- TwoSampleMR::add_metadata(out, cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")) %>%
    # dplyr::select(rsid=SNP, chr, position=pos, id=id.outcome, beta=beta.outcome, se=se.outcome, p=pval.exposure, ea=effect_allele.outcome, nea=other_allele.outcome, eaf=eaf.outcome, units=units.outcome, samplesize=contains("size")) %>%
    as.data.frame())
  self$instrument_outcome <- out %>% dplyr::arrange(., chr, SNP)
  invisible(self)
})


#' @description
#' This function harmonises the alleles and effects between the exposure and outcome.
#' @param exp Intsruments for the exposure that are selected by using the provided methods in CAMERA (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param out Intsruments for the outcome by using \code{make_outcome_data()}. Default is x$instrument_outcome.
#' @return Data frame in x$harmonised_dat_sem
CAMERA$set("public", "harmonised_dat", function(exp = self$instrument_raw, out = self$instrument_outcome) {
  dx <- dplyr::inner_join(
    subset(exp, id == self$exposure_ids[[1]]),
    subset(exp, id == self$exposure_ids[[2]]),
    by = "rsid"
  ) %>%
    dplyr::select(SNP = rsid, x1 = beta.x, x2 = beta.y, xse1 = se.x, xse2 = se.y, p1 = p.x, p2 = p.y)
  dy <- dplyr::inner_join(
    subset(out, id.outcome == self$outcome_ids[[1]]),
    subset(out, id.outcome == self$outcome_ids[[2]]),
    by = "SNP"
  ) %>%
    dplyr::select(SNP = SNP, y1 = beta.outcome.x, y2 = beta.outcome.y, yse1 = se.outcome.x, yse2 = se.outcome.y)
  dat <- dplyr::inner_join(dx, dy, by = "SNP")
  self$harmonised_dat_sem <- dat
})
