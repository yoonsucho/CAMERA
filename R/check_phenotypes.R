#' @description
#'  This function evaluates how competitable genetic associations in population 1 are with those in population 2. The function checks whether 1) the chosen IDs for the exposure or outcome can be used for the further steps and 2) the units for the genetic associations are comparable between populations.
#' @param ids ID for the exposure or the outcome. Default is x$exposure_ids.
#' @return Table of the result.
CAMERA$set("public", "check_phenotypes", function(ids = self$exposure_ids) {
  o <- lapply(ids, function(i) {
    tryCatch(
      {
        suppressMessages(exp <- unique(TwoSampleMR::extract_instruments(outcomes = i)))
        other_ids <- ids[!ids %in% i]

        o <- lapply(other_ids, function(j) {
          suppressMessages(out <- TwoSampleMR::extract_outcome_data(snps = exp$SNP, outcomes = j))
          suppressMessages(d <- TwoSampleMR::harmonise_data(exp, out))

          res <- suppressMessages(TwoSampleMR::mr(d, method = "mr_ivw")) %>%
            {
              tibble::tibble(Reference = i, Replication = j, nsnp = .$nsnp, agreement = .$b, se = .$se, pval = .$pval)
            }

          message(paste0("Instrument associations between ", i, " and ", j, " is ", round(res$agreement, 3), "; NSNP=", res$nsnp))

          suppressMessages(
            t <- d %>% data.frame() %>%
              TwoSampleMR::add_metadata(., cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")) %>%
              dplyr::summarise(unit_ref = .$units.exposure[1], unit_rep = .$units.outcome[1])
          )

          if (any(is.na(t))) {
            message("Unit information is missing: See vignettes")
            print(t)
          }

          if (!any(is.na(t)) & t$unit_ref != t$unit_rep) {
            message("Units for the beta are different across the populations: Run x$standardise_data()")
            print(t)
          }
          return(res)
        })
        return(o %>% dplyr::bind_rows())
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })
  print(o %>% dplyr::bind_rows())
})
