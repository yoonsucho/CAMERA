#' R6 class for CAMERA
#'
#' @description
#' A simple wrapper function.
#' Using a summary set, identify set of instruments for the traits, and peform SEM MR to test the association across the population.
#' @export

CAMERA <- R6::R6Class("CAMERA", list(
  output = list(),
  source = NULL,
  exposure_ids = NULL,
  outcome_ids = NULL,
  exposure_metadata = NULL,
  outcome_metadata = NULL,
  radius = NULL,
  pops = NULL,
  bfiles = NULL,
  plink = NULL,
  clump_pop = NULL,
  instrument_raw = NULL,
  instrument_regions = NULL,
  ld_matrices = NULL,
  susie_results = NULL,
  paintor_results = NULL,
  mscaviar_results = NULL,
  expected_replications = NULL,
  instrument_region_zscores = NULL,
  instrument_fema = NULL,
  instrument_fema_regions = NULL,
  instrument_maxz = NULL,
  instrument_susie = NULL,
  instrument_paintor = NULL,
  instrument_mscaviar = NULL,
  harmonised_data_check = NULL,
  standardised_instrument_raw = NULL,
  standardised_instrument_maxz = NULL,
  standardised_instrument_susie = NULL,
  standardised_instrument_paintor = NULL,
  standardised_instrument_mscaviar = NULL,
  standardised_outcome = NULL,
  instrument_specificity = NULL,
  instrument_specificity_summary = NULL,
  instrument_outcome = NULL,
  instrument_outcome_regions = NULL,
  harmonised_dat_sem = NULL,
  harmonised_dat = NULL,
  sem_result = NULL,
  pleiotropy_agreement = NULL,
  pleiotropy_outliers = NULL,
  pleiotropy_Q_outliers = NULL,
  summary = NULL,
  mrres = NULL,
  instrument_heterogeneity_per_variant = NULL,
  mrgxe_res = NULL,

  # for convenience can migrate the results from a previous CAMERA into this one
  #' @description
  #' Migrate the results from a previous CAMERA
  #' @param x R6 Environment created for CAMERA. Default = x
  import = function(x) {
    nom <- names(self)
    for (i in nom)
    {
      if (!i %in% c(".__enclos_env__", "clone") & is.null(self[[i]])) {
        self[[i]] <- x[[i]]
      }
    }
  },

  assign = function(...) {
    l <- list(...)
    lapply(names(l), \(n) {
      message("Assigning ", n)
      try(self[[n]] <- l[[n]])
    })
  },

  import_from_local = function(instrument_raw, instrument_outcome, instrument_regions, instrument_outcome_regions, exposure_ids, outcome_ids, pops, ...) {
    self[["instrument_raw"]] <- instrument_raw %>% generate_vid()
    self[["instrument_outcome"]] <- instrument_outcome  %>% generate_vid()
    self[["instrument_regions"]] <- lapply(instrument_regions, \(x) lapply(x, generate_vid))
    self[["instrument_outcome_regions"]] <- lapply(instrument_outcome_regions, \(x) lapply(x, generate_vid))
    self[["exposure_ids"]] <- exposure_ids
    self[["outcome_ids"]] <- outcome_ids
    self$source <- "Local"
    self$pops <- pops
    # Get instrument_outcome

    self$assign(...)
  },

  # Methods
  #' @description
  #' Create a new dataset and initialise an R interface
  #' @param exposure_ids Exposures IDs obtained from IEU GWAS database (https://gwas.mrcieu.ac.uk/) for each population
  #' @param outcome_ids Outcome IDs obtained from IEU GWAS database (https://gwas.mrcieu.ac.uk/) for each population
  #' @param pops Ancestry information for each population (i.e. AFR, AMR, EUR, EAS, SAS)
  #' @param bfiles Locations of LD reference files for each population (Download from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz)
  #' @param plink Location of executable plink (ver.1.90 is recommended)
  #' @param radius Genomic window size to extract SNPs
  #' @param clump_pop Reference population for clumping
  #' @param x Import data where available
  initialize = function(exposure_ids = NULL, outcome_ids = NULL, pops = NULL, bfiles = NULL, plink = NULL, radius = NULL, clump_pop = NULL, x = NULL) {
    if (!is.null(x)) {
      import(x)
    }
    self$exposure_ids <- exposure_ids
    self$outcome_ids <- outcome_ids
    self$pops <- pops
    self$bfiles <- bfiles
    self$plink <- plink
    self$radius <- radius
    self$clump_pop <- clump_pop
    if(!is.null(exposure_ids)) {
      self$source <- "OpenGWAS"
      self$get_metadata()
    }
  }
))
