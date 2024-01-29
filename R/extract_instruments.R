generate_vid <- function(d, ea = "ea", nea = "nea", eaf = "eaf", beta = "beta", rsid = "rsid", chr = "chr", position = "position") {
  toflip <- d[[ea]] > d[[nea]]
  d[[eaf]][toflip] <- 1 - d[[eaf]][toflip]
  d[[beta]][toflip] <- d[[beta]][toflip] * -1
  temp <- d[[nea]][toflip]
  d[[nea]][toflip] <- d[[ea]][toflip]
  d[[ea]][toflip] <- temp
  d[["rsido"]] <- d[[rsid]]
  d[[rsid]] <- paste0(d[[chr]], ":", d[[position]], "_", strtrim(d[[ea]], 5), "_", strtrim(d[[nea]], 5))
  d
}


#' @description
#'  This function searches for GWAS significant SNPs (P < 5E-8) for a specified set of the exposures. This method is equivalant to the instrumnet extraction method for Multivariable MR. Reference here: https://mrcieu.github.io/TwoSampleMR/reference/mv_extract_exposures.html.
#' @param exposure_ids ID for the exposure. Default is x$exposure_ids.
#' @return Data frame in x$instrument_raw
#' @importFrom ieugwasr variants_rsid
CAMERA$set("public", "extract_instruments", function(exposure_ids = self$exposure_ids, ...) {
  suppressMessages(instrument_raw <- TwoSampleMR::mv_extract_exposures(exposure_ids, ...))
  # Add chromosome and position
  suppressMessages(instrument_raw <- TwoSampleMR::add_metadata(instrument_raw, cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")))
  suppressMessages(instrument_raw <- ieugwasr::variants_rsid(unique(instrument_raw$SNP)) %>%
    dplyr::select(SNP = query, chr, position = pos) %>%
    dplyr::inner_join(., instrument_raw, by = "SNP") %>%
    dplyr::arrange(id.exposure, chr, position))
  # Arrange to be in order of exposure_ids
  # Rename columns
  instrument_raw <- lapply(self$exposure_ids, function(id) {
    subset(instrument_raw, id.exposure == id)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(rsid = SNP, chr, position, id = id.exposure, beta = beta.exposure, se = se.exposure, p = pval.exposure, ea = effect_allele.exposure, nea = other_allele.exposure, eaf = eaf.exposure, units = units.exposure, samplesize = contains("size")) %>%
    dplyr::mutate(method = "raw") %>%
    as.data.frame()
  # Check if top hits are significant in both populations
  t <- instrument_raw %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(sum(p < 5e-8))
  id <- list()
  id <- t$id[t$`sum(p < 5e-08)` < 1]
  if (length(id) > 0) {
    message(paste0("Caution: No SNPs reached genome-wide significance threshold for the trait in ", id))
  }
  self$instrument_raw <- generate_vid(instrument_raw)
  invisible(self)
})

# Here the idea is that pop1 and pop2 might share an instrument, but the tophit for pop1 is not the causal variant
# Hence, in pop1 it is in LD with the causal variant but not in pop2
#' @description
#' This function extract genomic regions around each instrument (e.g. 50kb) obtained from \code{x$extract_instruments()}.
#' @param radius Set a range of the region to search
#' @param instrument_raw A set of instruments obtained from \code{x$extract_instruments()}
#' @param exposure_ids ID for the exposure. Default is x$exposure_ids.
#' @return Data frame in x$instrument_regions
#' @importFrom ieugwasr associations
CAMERA$set("public", "extract_instrument_regions", function(radius = self$radius, instrument_raw = self$instrument_raw, exposure_ids = self$exposure_ids) {
  # return a list of lists e.g.
  # region1:
  # pop1:
  # data frame
  # pop2:
  # data frame
  # create list of regions in chr:pos format
  temp <- subset(instrument_raw, !duplicated(paste(chr, position)))
  regions <- paste0(temp$chr, ":", temp$position - radius, "-", temp$position + radius)

  # Lookup each region in each exposure
  self$instrument_regions <- lapply(regions, function(r) {
    tryCatch(
      {
        message(r)
        a <- ieugwasr::associations(r, self$exposure_ids) %>%
          dplyr::arrange(position) %>%
          generate_vid()
        message(nrow(a))
        a <- lapply(self$exposure_ids, function(i) {
          subset(a, id == i) %>%
          dplyr::filter(!duplicated(rsid))
        })

        # subset to keep only the same SNPs across datasets
        # Make sure they have the same effect and other alleles
        rsids <- Reduce(intersect, lapply(a, function(x) x$rsid))
        a <- lapply(a, function(x) {
          subset(x, rsid %in% rsids)
        })

        # check effect and other alleles are the same
        ea <- a[[1]]$ea
        a <- lapply(a, function(x) {
          index <- x$ea != ea
          if (sum(index) > 0) {
            x$beta[index] <- x$beta[index] * -1
            nea <- x$nea[index]
            x$nea[index] <- x$ea[index]
            x$ea[index] <- x$nea[index]
            x$eaf[index] <- 1 - x$eaf[index]
            x <- subset(x, nea == a[[1]]$nea)
          }
          return(x)
        })
        rsids <- Reduce(intersect, lapply(a, function(x) x$rsid))
        a <- lapply(a, function(x) {
          subset(x, rsid %in% rsids) %>%
            dplyr::arrange(chr, position)
        })
        names(a) <- self$exposure_ids
        return(a)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })
  names(self$instrument_regions) <- temp$rsid
})


#' @description
#' The function searches for a SNP that is best associated across the poupulations within the identified regions by using \code{x$extract_instrument_regions()}.
#' @param instrument_raw Instruments for the exposures by using \code{x$extract_instrument()}
#' @param instrument_regions Genomic regions identified by using \code{x$extract_instrument_regions()}
#' @return List of z scores for the chosen SNPs in x$instrument_region_zscores. Data frame in x$instrument_maxz
CAMERA$set("public", "scan_regional_instruments", function(instrument_raw = self$instrument_raw, instrument_regions = self$instrument_regions) {
  # Simple method to choose the best SNP in the region by
  # - normalising the Z scores for each trait (to be in range -1 to 1)
  # - adding the z scores together across traits
  # - choosing the largest abs(z) across all traits
  instrument_regions <- instrument_regions[lengths(instrument_regions) != 0]
  temp <- lapply(instrument_regions, function(r) {
    lapply(r, function(id) {
      id$beta / id$se
    }) %>% dplyr::bind_cols()
  })
  z <- lapply(temp, function(x) {
    apply(x, 2, function(y) y / max(abs(y))) %>%
      rowSums()
  })
  self$instrument_region_zscores <- lapply(1:length(temp), function(i) {
    temp[[i]] %>%
      dplyr::mutate(
        zsum = z[[i]],
        rsid = instrument_regions[[i]][[1]]$rsid,
        chr = instrument_regions[[i]][[1]]$chr,
        position = instrument_regions[[i]][[1]]$position,
      )
  })
  names(self$instrument_region_zscores) <- names(instrument_regions)
  self$instrument_maxz <- names(self$instrument_region_zscores) %>%
    lapply(., function(r) {
      o <- self$instrument_region_zscores[[r]] %>%
        dplyr::arrange(desc(abs(zsum))) %>%
        dplyr::slice(1) %>%
        dplyr::mutate(original = r)
      instrument_regions[[r]] %>%
        lapply(., function(id) {
          subset(id, rsid == o$rsid) %>%
            dplyr::mutate(original_rsid = r, zsum = o$zsum)
        }) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(id, chr, position)
  t <- instrument_raw %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::select(id, units, samplesize)
  self$instrument_maxz <- dplyr::left_join(self$instrument_maxz, t, by = "id") %>% as.data.frame()
  self$instrument_maxz <- lapply(self$exposure_ids, function(i) {
    subset(self$instrument_maxz, id == i) %>% dplyr::distinct(., rsid, .keep_all = TRUE)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(method = "maxz")
})


#' @description
#' The function draws a plot shows how CAMERA selected the instruments for the exposure from the selected genomic regions.
#' @param region X axis. Number of genomic regions to be shown in the plot. Default is 10.
#' @param instrument_region_zscores Y axis. Z scores based on the SNP-exposure associations.
#' @param instruments Use this option to draw a separate plot for the selcted instruments.
#' @param comparison Use this option to compare the selected instruments by different instrument selection methods in one plot.
#' @return Plot
#' @importFrom ggplot2 ggplot aes geom_point facet_grid geom_smooth scale_colour_brewer scale_x_log10 scale_y_log10 xlab ylab
CAMERA$set("public", "plot_regional_instruments_maxz", function(instrument_region_zscores = self$instrument_region_zscores, instruments = self$instrument_raw, region = 1:min(10, nrow(instruments)), comparison = FALSE) {
  a <- instrument_region_zscores[region]
  a <- names(a) %>%
    lapply(., function(n) {
      o <- dplyr::bind_rows(a[[n]])
      o$original_rsid <- n
      return(o)
    }) %>%
    dplyr::bind_rows()

  colnum <- which(names(a) == "zsum")

  nom <- names(a)[1:colnum]

  a <- tidyr::pivot_longer(a, nom)

  if (comparison == FALSE) {
    a$selected <- a$rsid %in% instruments$rsid
    p <- ggplot2::ggplot(a, ggplot2::aes(x = position, y = value)) +
      ggplot2::geom_point(ggplot2::aes(colour = name)) +
      ggplot2::geom_point(data = subset(a, selected), colour = "black") +
      ggplot2::facet_grid(name ~ original_rsid, scale = "free")
  }

  if (comparison == TRUE) {
    a$selected_raw <- a$rsid %in% self$instrument_raw$rsid
    a$selected_maxz <- a$rsid %in% self$instrument_maxz$rsid
    a$selected_susie <- a$rsid %in% self$instrument_susie$rsid
    a$selected_paintor <- a$rsid %in% self$instrument_paintor$rsid
    p <- ggplot2::ggplot(a, ggplot2::aes(x = position, y = value)) +
      ggplot2::geom_point(ggplot2::aes(colour = name)) +
      ggplot2::geom_point(data = subset(a, selected_raw), colour = "black") +
      ggplot2::geom_point(data = subset(a, selected_maxz), colour = "purple", shape = 15) +
      ggplot2::geom_point(data = subset(a, selected_susie), colour = "orange", shape = 17) +
      ggplot2::geom_point(data = subset(a, selected_paintor), colour = "brown", shape = 18) +
      ggplot2::facet_grid(name ~ original_rsid, scale = "free")
  }
  p
})

#' Generate LD matrices for instrument regions
#'
#' @description
#' If we want to do fine mapping we need to get an LD matrix for the whole region (for each population)
#' We then need to harmonise the LD matrix to the summary data, and the summary datasets to each other
#' The fuction obtains an LD matrix for the selected genomic regions.
#'
#' @param instrument_regions Genomic regions identified by using \code{x$extract_instrument_regions()}
#' @param bfiles Location of LD reference files for each population (Download from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz)
#' @param pops Ancestry information for each population (i.e. AFR, AMR, EUR, EAS, SAS)
#' @param plink Location of executable plink (ver.1.90 is recommended)
#' @return Data frame of LD matrix (x$ld_matrices)
#' @importFrom ieugwasr ld_matrix
CAMERA$set("public", "regional_ld_matrices", function(instrument_regions = self$instrument_regions, bfiles = self$bfiles, pops = self$pops, plink = self$plink) {
  if (!is.null(bfiles)) {
    stopifnot(length(bfiles) == length(self$exposure_ids))
  }
  if (!is.null(pops)) {
    stopifnot(length(pops) == length(self$exposure_ids))
  }

  regions <- names(instrument_regions)
  ld_matrices <- lapply(regions, function(r) {
    tryCatch(
      {
        d <- instrument_regions[[r]]
        exp <- self$exposure_ids
        o <- lapply(1:length(self$exposure_ids), function(i) {
          ld <- ieugwasr::ld_matrix(d[[exp[i]]]$rsid, pop = pops[i], bfile = bfiles[i], plink = self$plink, with_alleles = TRUE)
          code1 <- paste0(d[[exp[i]]]$rsid, "_", d[[exp[i]]]$ea, "_", d[[exp[i]]]$nea)
          code2 <- paste0(d[[exp[i]]]$rsid, "_", d[[exp[i]]]$nea, "_", d[[exp[i]]]$ea)
          rem_index <- !(code1 %in% colnames(ld) | code2 %in% colnames(ld))
          flip_index <- !colnames(ld) %in% code1
          if (any(rem_index)) {
            rem <- d[[exp[i]]]$rsid[rem_index]
          }
          if (!any(rem_index)) {
            rem <- NULL
          }
          if (any(flip_index)) {
            message("Flipping ", sum(flip_index))
            ld[flip_index, ] <- ld[flip_index, ] * -1
            ld[, flip_index] <- ld[, flip_index] * -1
          }
          return(list(ld = ld, rem = rem))
        })

        rem <- lapply(o, function(x) x$rem) %>%
          unlist() %>%
          unique()
        o <- lapply(o, function(x) x$ld)
        names(o) <- exp
        for (i in exp)
        {
          instrument_regions[[1]][[i]]$ld_unavailable <- instrument_regions[[1]][[i]]$rsid %in% rem
        }
        o <- lapply(o, function(ld) {
          rs <- strsplit(colnames(ld), "_") %>% sapply(., function(x) x[1])
          ld[!rs %in% rem, !rs %in% rem]
        })
        return(o)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })
  names(ld_matrices) <- regions
  self$ld_matrices <- ld_matrices
  invisible(self)
})


CAMERA$set("private", "greedy_remove", function(r, thresh) {
  diag(r) <- 0
  r <- abs(r)
  flag <- 1
  rem <- c()
  nom <- colnames(r)
  while (flag == 1) {
    count <- apply(r, 2, function(x) sum(x >= thresh))
    if (any(count > 0)) {
      worst <- which.max(count)[1]
      rem <- c(rem, names(worst))
      r <- r[-worst, -worst]
    } else {
      flag <- 0
    }
  }
  return(which(nom %in% rem))
})
