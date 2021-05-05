#' R6 class for MR.TRAM
#'
#' @description
#' A simple wrapper function.
#' Using a summary set, identify set of instruments for the traits, and peform SEM MR to test the association across the population.
#' @export

MultiAncestrySummarySet <- R6::R6Class("MultiAncestrySummarySet", list(
  output = list(),
  exposure_ids=NULL,
  outcome_ids=NULL,
  radius=NULL,
  pops=NULL,
  bfiles=NULL,
  plink=NULL,
  clump_pop=NULL,
  instrument_raw=NULL,
  instrument_regions = NULL,
  ld_matrices = NULL,
  susie_results = NULL,
  paintor_results = NULL,
  expected_replications = NULL,
  instrument_region_zscores = NULL,
  instrument_maxz = NULL,
  instrument_specificity = NULL,
  instrument_specificity_summary = NULL,
  instrument_outcome = NULL,
  harmonised_dat = NULL,
  sem_result = NULL,

# for convenience can migrate the results from a previous MultiAncestrySummarySet into this one
#' @description
#' Migrate the results from a previous MultiAncestrySummarySet
  import = function(x) {
    nom <- names(self)
    for(i in nom)
    {
      if(! i %in% c(".__enclos_env__", "clone") & is.null(self[[i]]))
      {
        self[[i]] <- x[[i]]
      }
    }
  },

# Methods
#' @description
#' Create a new dataset and initialise an R interface
#' @param exposure_ids IDs for the exposures from each population.
#' @param outcome_ids IDs for the outcomes from each population.
#' @param pops Ancestry information; AFR, AMR, EUR, EAS, SAS 
#' @param bfiles Location of LD reference file (Download from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz)
#' @param plink Location of executable plink (ver.1.90)
#' @param radius Set a range of the region to look for
#' @param clump_pop insert
#' @param x insert
  initialize = function(exposure_ids=NULL, outcome_ids=NULL, pops=NULL, bfiles=NULL, plink=NULL, radius=NULL, clump_pop=NULL, x=NULL)
  {
    if(!is.null(x))
    {
      import(x)
    }
    self$exposure_ids <- exposure_ids
    self$outcome_ids <- outcome_ids
    self$pops <- pops
    self$bfiles <- bfiles
    self$plink <- plink
    self$radius <- radius
    self$clump_pop <- clump_pop
  },

  # e.g. could just use mv_extract_exposures
  # - for each exposure_id it identifies the instruments
  # - then it tries to make these independent
  # - then looks up the same set of instruments in all exposures
#' @description
#'  Identifies the instruments for the exposure
#' @param exposure_ids ID for the exposure. Default is x$exposure_ids.

  extract_instruments = function(exposure_ids=self$exposure_ids, ...)
  {
    # Use MVMR method to do the initial extraction
    # It gets the tophits in each exposure
    # Then randomly chooses one SNP per region to keep using clumping
    instrument_raw <- TwoSampleMR::mv_extract_exposures(exposure_ids, ...)

    # Add chromosome and position
    instrument_raw <- ieugwasr::variants_rsid(unique(instrument_raw$SNP)) %>%
      dplyr::select(SNP=query, chr, position=pos) %>%
      dplyr::inner_join(., instrument_raw, by="SNP") %>%
      dplyr::arrange(id.exposure, chr, position)

    # Arrange to be in order of exposure_ids
    # Rename columns
    instrument_raw <- lapply(self$exposure_ids, function(id) {
      subset(instrument_raw, id.exposure==id)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::select(rsid=SNP, chr, position, id=id.exposure, beta=beta.exposure, se=se.exposure, p=pval.exposure, ea=effect_allele.exposure, nea=other_allele.exposure, eaf=eaf.exposure)
    instrument_raw %>% as.data.frame
    self$instrument_raw <- instrument_raw
    invisible(self)
  },

  # Here the idea is that pop1 and pop2 might share an instrument, but the tophit for pop1 is not the causal variant
  # Hence, in pop1 it is in LD with the causal variant but not in pop2
  # So we extract a region around each instrument (e.g. 50kb)
  # Search for a SNP that is best associated in both pop1 and pop2
#' @description
#' Extract a region around each instrument obtained from \code{x$extract_instruments()}.
#' @param radius Set a range of the region to search
#' @param instrument_raw A set of instruments obtained from \code{x$extract_instruments()}
#' @param exposure_ids ID for the exposure. Default is x$exposure_ids.
  extract_instrument_regions = function(radius=self$radius, instrument_raw=self$instrument_raw, exposure_ids=self$exposure_ids)
  {
    # return a list of lists e.g.
    # region1:
    # pop1:
    # data frame
    # pop2:
    # data frame

    # create list of regions in chr:pos format
    regions <- unique(paste0(instrument_raw$chr, ":", instrument_raw$position-radius, "-", instrument_raw$position+radius))

    # Lookup each region in each exposure
    self$instrument_regions <- lapply(regions, function(r) {
      message(r)
      a <- lapply(self$exposure_ids, function(id) {
        message(id)
        ieugwasr::associations(r, id) %>%
          dplyr::arrange(position) %>%
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
        if(sum(index) > 0)
        {
          x$beta[index] <- x$beta[index] * -1
          nea <- x$nea[index]
          x$nea[index] <- x$ea[index]
          x$ea[index] <- x$nea[index]
          x$eaf[index] <- 1-x$eaf[index]
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
    })
    names(self$instrument_regions) <- unique(instrument_raw$rsid)
  },


  # Once a set of instruments is chosen we can ask
  # - what fraction of those primarily identified in pop1 replicate in pop2?
  # - what fraction of those primarily identified in pop1 have the same sign as in pop2?
  # - compare these to what is expected by chance under the hypothesis that the effect estimates are the same
  # this will provide some evidence for whether lack of replication is due to (e.g.) GxE
#' @description
#' Estimate what fraction of the instuments is expected to be replicated under the hypothesis that the effect estimates are the same
#' @param instrument A set of instruments obtained from \code{x$extract_instruments()} or \code{scan_regional_instruments()}
#' @param alpha Significance level. Default is 0.05/number of SNPs (Bonferroni)

  estimate_instrument_specificity = function(instrument, alpha="bonferroni")
  {
    if(alpha=="bonferroni")
    {
      alpha <- 0.05/nrow(instrument)
    }
    o <- lapply(self$exposure_ids, function(i)
    {
      m <- subset(instrument, id == i & p < 5e-8)
      other_ids <- self$exposure_ids[!self$exposure_ids %in% i]
      o <- lapply(other_ids, function(j)
      {
        n <- subset(instrument, id == j & rsid %in% m$rsid)
        o <- n %>%
          {self$prop_overlap(m$beta, .$beta, m$se, .$se, alpha)}
        o$res <- o$res %>%
          dplyr::mutate(discovery=i, replication=j) %>%
          dplyr::select(discovery, replication, dplyr::everything())
        o$variants <- o$variants %>%
          dplyr::mutate(
            discovery=i,
            replication=j,
            rsid=n$rsid,
            p=n$p,
            distinct=sig > 0.8 & p > 0.1
          ) %>%
          dplyr::select(discovery, replication, dplyr::everything())
        return(o)
      })
      overall <- lapply(o, function(x) { x$res }) %>% dplyr::bind_rows()
      pervariant <- lapply(o, function(x) { x$variants }) %>% dplyr::bind_rows()
      return(list(overall=overall, pervariant=pervariant))
    })
    self$instrument_specificity_summary <- lapply(o, function(x) x$overall) %>% dplyr::bind_rows()
    self$instrument_specificity <- lapply(o, function(x) x$pervariant) %>% dplyr::bind_rows()
    return(self$instrument_specificity_summary)
  },

#' @description
#' insert
#' @param b_disc insert
#' @param b_rep insert
#' @param se_disc insert
#' @param se_rep insert
#' @param alpha insert
  prop_overlap = function(b_disc, b_rep, se_disc, se_rep, alpha)
  {
    p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
    p_sig <- pnorm(-abs(b_disc) / se_disc + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_disc - qnorm(alpha / 2)))
    p_rep <- pnorm(abs(b_rep)/se_rep, lower.tail=FALSE)
    res <- dplyr::tibble(
      nsnp=length(b_disc),
      metric=c("Sign", "Sign", "P-value", "P-value"),
      datum=c("Expected", "Observed", "Expected", "Observed"),
      value=c(sum(p_sign, na.rm=TRUE), sum(sign(b_disc) == sign(b_rep)), sum(p_sig, na.rm=TRUE), sum(p_rep < alpha, na.rm=TRUE))
    )
    return(list(res=res, variants=dplyr::tibble(sig=p_sig, sign=p_sign)))
  },


#' @description
#' insert
#' @param instrument_raw insert
#' @param instrument_regions insert
  scan_regional_instruments = function(instrument_raw=self$instrument_raw, instrument_regions=self$instrument_regions)
  {
    # Simple method to choose the best SNP in the region by
    # - normalising the Z scores for each trait (to be in range -1 to 1)
    # - adding the z scores together across traits
    # - choosing the largest abs(z) across all traits

    temp <- lapply(self$instrument_regions, function(r) {
      lapply(r, function(id) {
        id$beta / id$se
      }) %>% dplyr::bind_cols()
    })
    z <- lapply(temp, function(x)
    {
      apply(x, 2, function(y) y/max(abs(y))) %>%
        rowSums()
    })
    self$instrument_region_zscores <- lapply(1:length(temp), function(i)
    {
      temp[[i]] %>%
        dplyr::mutate(
          zsum = z[[i]],
          rsid = self$instrument_regions[[i]][[1]]$rsid,
          chr = self$instrument_regions[[i]][[1]]$chr,
          position = self$instrument_regions[[i]][[1]]$position,
        )
    })
    names(self$instrument_region_zscores) <- names(self$instrument_regions)
    self$instrument_maxz <- names(self$instrument_region_zscores) %>%
      lapply(., function(r) {
        o <- self$instrument_region_zscores[[r]] %>%
          dplyr::arrange(desc(abs(zsum))) %>%
          dplyr::slice(1) %>%
          dplyr::mutate(original=r)
        self$instrument_regions[[r]] %>%
          lapply(., function(id){
            subset(id, rsid == o$rsid) %>%
              dplyr::mutate(original_rsid = r, zsum=o$zsum)
          }) %>% dplyr::bind_rows()
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(id, chr, position)
    self$instrument_maxz <- lapply(self$exposure_ids, function(i)
    {
      subset(self$instrument_maxz, id == i)
    }) %>% dplyr::bind_rows()

  },

#' @description
#' insert
#' @param region insert
#' @param instrument_region_zscores insert
#' @param instruments insert
  plot_regional_instruments = function(region=1:min(10, nrow(instruments)), instrument_region_zscores=self$instrument_region_zscores, instruments=self$instrument_raw)
  {
    a <- instrument_region_zscores[region]
    a <- names(a) %>% lapply(., function(n)
    {
      o <- dplyr::bind_rows(a[[n]])
      o$original_rsid <- n
      return(o)
    }) %>% dplyr::bind_rows()
    colnum <- which(names(a) == "zsum")
    nom <- names(a)[1:colnum]
    a <- tidyr::pivot_longer(a, nom)
    a$selected <- a$rsid %in% instruments$rsid
    ggplot2::ggplot(a, ggplot2::aes(x=position, y=value)) +
      ggplot2::geom_point(ggplot2::aes(colour=name)) +
      ggplot2::geom_point(data=subset(a, selected), , colour="black") +
      ggplot2::facet_grid(name ~ original_rsid, scale="free")
  },


  # to build on extract_instrument_regions we can do finemapping
  # If we want to do fine mapping we need to get an LD matrix for the whole region (for each population)
  # We then need to harmonise the LD matrix to the summary data, and the summary datasets to each other
#' @description
#' insert
#' @param instrument_regions insert
#' @param bfiles insert
#' @param pops insert
#' @param plink insert
  regional_ld_matrices = function(instrument_regions=self$instrument_regions, bfiles=self$bfiles, pops=self$pops, plink=self$plink)
  {

    if(!is.null(bfiles))
    {
      stopifnot(length(bfiles) == length(self$exposure_ids))
    }
    if(!is.null(pops))
    {
      stopifnot(length(pops) == length(self$exposure_ids))
    }

    regions <- names(instrument_regions)
    ld_matrices <- lapply(regions, function(r)
    {
      d <- instrument_regions[[1]]
      o <- lapply(1:length(self$exposure_ids), function(i)
      {
        ld <- ieugwasr::ld_matrix(d[[self$exposure_ids[1]]]$rsid, pop=self$pops[1], bfile=self$bfiles[1], plink=self$plink, with_alleles=TRUE)
        code1 <- paste0(d[[self$exposure_ids[i]]]$rsid, "_", d[[self$exposure_ids[i]]]$ea, "_", d[[self$exposure_ids[i]]]$nea)
        code2 <- paste0(d[[self$exposure_ids[i]]]$rsid, "_", d[[self$exposure_ids[i]]]$ea, "_", d[[self$exposure_ids[i]]]$nea)
        rem_index <- ! (code1 %in% colnames(ld) | code2 %in% colnames(ld))
        flip_index <- ! colnames(ld) %in% code1
        if(any(rem_index))
        {
          rem <- d[[self$exposure_ids[i]]]$rsid[rem_index]
        }
        if(any(flip_index))
        {
          message("Flipping ", sum(flip_index))
          ld[flip_index, ] <- ld[flip_index, ] * -1
          ld[, flip_index] <- ld[, flip_index] * -1
        }
        return(list(ld=ld, rem=rem))
      })
      rem <- lapply(o, function(x) x$rem) %>% unlist() %>% unique()
      o <- lapply(o, function(x) x$ld)
      names(o) <- self$exposure_ids
      for(i in self$exposure_ids)
      {
        self$instrument_regions[[r]][[i]]$ld_unavailable <- self$instrument_regions[[r]][[i]]$rsid %in% rem
      }
      o <- lapply(o, function(ld)
      {
        rs <- strsplit(colnames(ld), "_") %>% sapply(., function(x) x[1])
        ld[! rs %in% rem, ! rs %in% rem]
      })
      return(o)
    })
    names(ld_matrices) <- regions
    self$ld_matrices <- ld_matrices
    invisible(self)
  },


  # for each instrument region + ld matrix we can perform susie finemapping
  # do this independently in each population
  # find credible sets that overlap - to try to determine the best SNP in a region to be used as instrument
  susie_finemap_regions = function() {},

  # PAINTOR allows finemapping jointly across multiple populations
  # returns a posterior probability of inclusion for each SNP
  # Choose the variant with highest posterior probability and associations in each exposure
  paintor_finemap_regions = function() {},


  # Once a set of instruments is defined then extract the outcome data
  # Could use
  # - raw results from extract_instruments
  # - maximised associations from extract_instrument_regions
  # - finemapped hits from susie_finemap_regions
  # - finemapped hits from paintor_finemap_regions
#' @description
#' insert
#' @param exp insert
#' @param p_exp insert
 extract_outcome_data = function(exp=self$instrument_raw, p_exp=1) {
    exp <- exp %>% as.data.frame
    dx <- dplyr::inner_join(
      subset(exp, id == self$exposure_ids[[1]]),
      subset(exp, id == self$exposure_ids[[2]]),
      by="rsid"
    ) %>%
      dplyr::select(SNP=rsid, x1=beta.x, x2=beta.y, xse1=se.x, xse2=se.y, p1=p.x, p2=p.y) %>%
      dplyr::filter(p1 < p_exp & p2 < p_exp)
    out <- TwoSampleMR::extract_outcome_data(snps=dx$SNP, outcomes=self$outcome_ids)

    dy <- dplyr::inner_join(
      subset(out, id.outcome == self$outcome_ids[[1]]),
      subset(out, id.outcome == self$outcome_ids[[2]]),
      by="SNP"
      ) %>%
        dplyr::select(SNP=SNP, y1=beta.outcome.x, y2=beta.outcome.y, yse1=se.outcome.x, yse2=se.outcome.y)
    dat <- dplyr::inner_join(dx, dy, by="SNP")   

    self$instrument_outcome <- out
    self$harmonised_dat <- dat

    invisible(self)
  },

  # Generate harmonised dataset

  # Perform basic SEM analysis of the joint estimates in multiple ancestries
#' @description
#' insert
#' @param harmonised_dat insert
  perform_basic_sem = function(harmonised_dat = self$harmonised_dat) {
    d <- harmonised_dat %>%
         dplyr::mutate(r1 = y1/x1) %>%
         dplyr::mutate(r2 = y2/x2) %>%
         dplyr::mutate(w1 = sqrt(x1^2 / yse1^2)) %>%
         dplyr::mutate(w2 = sqrt(x2^2 / yse2^2)) %>%
         dplyr::mutate(o1 = r1 * w1) %>%
         dplyr::mutate(o2 = r2 * w2)

    out <- list()
    out$ivw1 <- TwoSampleMR::mr_ivw(d$x1, d$y1, d$xse1, d$yse1) %>%
                  {tibble::tibble(model="IVW", pop=1, bivhat=.$b, se=.$se, pval=.$pval)}
    out$ivw2 <- TwoSampleMR::mr_ivw(d$x2, d$y2, d$xse2, d$yse2) %>%
                  {tibble::tibble(model="IVW", pop=2, bivhat=.$b, se=.$se, pval=.$pval)}
    out$rm1 <- summary(lm(o1 ~ -1 + w1, data=d)) %>%
                  {tibble::tibble(model="RadialIVW", pop=1, bivhat=.$coef[1,1], se=.$coef[1,2], pval=.$coef[1,4])}
    out$rm2 <- summary(lm(o2 ~ -1 + w2, data=d)) %>%
                  {tibble::tibble(model="RadialIVW", pop=2, bivhat=.$coef[1,1], se=.$coef[1,2], pval=.$coef[1,4])}

    runsem <- function(model, data, modname)
        {
          mod <- lavaan::sem(model, data=data)
          invisible(capture.output(mod <- lavaan::summary(mod, fit.measures=TRUE)))
          tibble::tibble(
                model=modname,
                pop=1:2,
                bivhat=mod$PE$est[1:2],
                se=mod$PE$se[1:2],
                pval=mod$PE$pvalue[1:2],
                aic=mod$FIT['aic']
              )
        }

    out$semA <- runsem('
                       y1 ~ biv*x1
                       y2 ~ biv*x2
                       ', d, "UnweightedSEMa")
    out$semB <- runsem('
                       y1 ~ biv_1*x1
                       y2 ~ biv_2*x2
                       ', d, "UnweightedSEMb")
    
    out$modC <- runsem('
                       o1 ~ biv*w1
                       o2 ~ biv*w2
                       ', d, "RadialSEMa")
    
    out$modD <- runsem('
                       o1 ~ biv_1*w1
                       o2 ~ biv_2*w2
                       ', d, "RadialSEMb")

    invisible(self$sem_result <- out %>% dplyr::bind_rows())
    print(self$sem_result)
  }

  # Plots
  # Plot regional associations for each instrument and each population
  # Plot finemapping results against regional association data

))









susie_overlaps <- function(su1, su2)
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
        l[[k]] <- tibble(s1=i, s2=j, v = v)
        k <- k + 1
      }
    }
  }
  l <- bind_rows(l)
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
}


greedy_remove <- function(r, thresh)
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
}

susie_finemap_region <- function(chr, position, radius, id1, id2, plink=NULL, bfile1=NULL, bfile2=NULL, pop1=NULL, pop2=NULL)
{
  r <- paste0(chr, ":", max(0, position-radius), "-", position+radius)
  message("Analysing ", r)
  message("Extracting summary data and LD for ", id1)
  a1 <- ieugwasr::associations(r, id1) %>%
    dplyr::arrange(position) %>%
    dplyr::filter(!duplicated(rsid)) %>%
    gwasglue::ieugwasr_to_TwoSampleMR(.)
  a1 <- TwoSampleMR::harmonise_ld_dat(a1, suppressMessages(suppressWarnings(ieugwasr::ld_matrix(a1$SNP, pop=pop1, plink=plink, bfile=bfile1, with_alleles=TRUE))))
  message("Extracting summary data and LD for ", id2)
  a2 <- ieugwasr::associations(r, id2) %>%
    dplyr::arrange(position) %>%
    dplyr::filter(!duplicated(rsid)) %>%
    gwasglue::ieugwasr_to_TwoSampleMR(.)
  a2 <- TwoSampleMR::harmonise_ld_dat(a2, suppressMessages(suppressWarnings(ieugwasr::ld_matrix(a2$SNP, pop=pop2, plink=plink, bfile=bfile2, with_alleles=TRUE))))

  snps <- a1$x$SNP[a1$x$SNP %in% a2$x$SNP]
  message(length(snps), " variants in common")
  index1 <- which(a1$x$SNP %in% snps)
  index2 <- which(a2$x$SNP %in% snps)
  a1$x <- a1$x[index1,]
  a1$ld <- a1$ld[index1, index1]
  a2$x <- a2$x[index2,]
  a2$ld <- a2$ld[index2, index2]

  stopifnot(all(a1$x$SNP == a2$x$SNP))

  message("Running susie for ", id1)
  su1 <- susieR::susie_rss(z=a1$x$beta.exposure / a1$x$se.exposure, R=a1$ld, check_R=FALSE)
  message("Running susie for ", id2)
  su2 <- susieR::susie_rss(z=a2$x$beta.exposure / a2$x$se.exposure, R=a2$ld, check_R=FALSE)

  message("Finding overlaps")
  suo <- susie_overlaps(su1, su2)
  out <- list(
    chr=chr,
    position=position,
    radius=radius,
    a1=a1$x,
    a2=a2$x,
    su1=su1,
    su2=su2,
    suo=suo
  )

  null <- c(is.null(su1$sets$cs), is.null(su2$sets$cs))
  if(null[1] & !null[2])
  {
    out$type <- "pop2"
  } else if(null[2] & !null[1]) {
    out$type <- "pop1"
  } else if(!null[1] & !null[2]) {
    out$type <- "shared"
  } else {
    out$type <- "drop"
  }

  if(out$type == "shared")
  {
    if(nrow(suo) == 0)
    {
      out$cs_overlap <- FALSE
      temp <- dplyr::inner_join(a1$x, a2$x, by="SNP") %>%
        dplyr::mutate(pvalrank = rank(pval.exposure.x) / nrow(a1$x) + rank(pval.exposure.y) / nrow(a1$x)) %>%
        dplyr::arrange(pvalrank)
      out$bestsnp <- temp$SNP[1]
    } else {
      out$cs_overlap <- TRUE
      out$bestsnp <- a1$x$SNP[suo$v[1]]
    }
  }

  if(out$type == "pop1")
  {
    out$cs_overlap <- FALSE
    out$bestsnp <- a1$x %>% dplyr::arrange(pval.exposure) %>% {.$SNP[1]}
  }
  if(out$type == "pop2")
  {
    out$cs_overlap <- FALSE
    out$bestsnp <- a2$x %>% dplyr::arrange(pval.exposure) %>% {.$SNP[1]}
  }
  return(out)
}

plot_region <- function(region)
{
  lapply(list(region$a1, region$a2), function(x)
    x %>%
      dplyr::select(pval.exposure, position.exposure, id.exposure)
  ) %>% dplyr::bind_rows() %>%
    ggplot2::ggplot(., ggplot2::aes(x=position.exposure, y=-log10(pval.exposure))) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(id.exposure ~ .)
}

#' Analyse regional data with PAINTOR
#'
#' @param regiondata Output from extract_regional_data
#' @param PAINTOR Path to PAINTOR executable. Default="PAINTOR"
#' @param workdir Working directory. Default=tempdir()
#'
#' @export
#' @return Results table with posterior inclusion probabilities
run_PAINTOR <- function(regiondata, PAINTOR="PAINTOR", workdir=tempdir())
{
  nid <- length(regiondata)
  snps <- regiondata[[1]]$x$SNP
  zs <- lapply(1:nid, function(i)
  {
    o <- list()
    o[[paste0("ZSCORE.P", i)]] <- regiondata[[i]]$x$beta.exposure / regiondata[[i]]$x$se.exposure
    return(as_tibble(o))
  }) %>% bind_cols()

  locus <- tibble(CHR = regiondata[[1]]$x$chr.exposure, POS = regiondata[[1]]$x$position.exposure, RSID = regiondata[[1]]$x$SNP) %>%
    bind_cols(., zs)
  write.table(locus, file=file.path(workdir, "Locus1"), row=F, col=T, qu=F)
  lapply(1:nid, function(i)
  {
    write.table(regiondata[[i]]$ld, file=file.path(workdir, paste0("Locus1.LD", i)), row=FALSE, col=FALSE, qu=FALSE)
  })
  write.table(tibble(null=rep(1, length(snps))), file=file.path(workdir, "Locus1.annotations"), row=F, col=T, qu=F)

  write.table(c("Locus1"), file=file.path(workdir, "input.files"), row=F, col=F, qu=F)

  LDname <- paste(paste0("LD", 1:nid), collapse=",")
  Zhead <- paste(names(zs), collapse=',')
  glue("{PAINTOR} -input {file.path(workdir,'input.files')} -Zhead {Zhead} -LDname {LDname} -in {workdir} -out {workdir} -mcmc -annotations null") %>%
    system()
  res <- data.table::fread(file.path(workdir, "Locus1.results"))
  unlink(workdir)
  return(res)
}


#' Plot PAINTOR results
#'
#' @param res Output from run_PAINTOR
#'
#' @export
#' @return plot
plot_PAINTOR <- function(res)
{
  res %>% pivot_longer(cols=c(ZSCORE.P1, ZSCORE.P2, Posterior_Prob)) %>%
    ggplot(., aes(x=POS, y=value)) +
    geom_point(aes(colour=name)) +
    facet_grid(name ~ ., scale="free_y")
}

#' Analyse regional data with MsCAVIAR
#'
#' @param regiondata Output from extract_regional_data
#' @param MsCAVIAR Path to MsCAVIAR executable. Default="MsCAVIAR"
#' @param workdir Working directory. Default=tempdir()
#'
#' @export
#' @return Results table with posterior inclusion probabilities
run_MsCAVIAR <- function(regiondata, n, MsCAVIAR="MsCAVIAR", workdir=tempdir())
{
  nid <- length(regiondata)
  snps <- regiondata[[1]]$x$SNP
  lapply(1:nid, function(i)
  {
    regiondata[[i]]$x %>%
      dplyr::mutate(z=beta.exposure/se.exposure) %>%
      dplyr::select(SNP, z) %>%
      write.table(., file=file.path(workdir, paste0("z.", i)), row=F, col=F, qu=F)
  })

  lapply(1:nid, function(i)
  {
    write.table(regiondata[[i]]$ld, file=file.path(workdir, paste0("ld.", i)), row=FALSE, col=FALSE, qu=FALSE)
  })

  thisdir <- getwd()
  setwd(workdir)
  write.table(paste0("ld.", 1:nid), file="ldfiles.txt", row=F, col=F, qu=F)
  write.table(paste0("z.", 1:nid), file="zfiles.txt", row=F, col=F, qu=F)
  glue("{MsCAVIAR} -l ldfiles.txt -z zfiles.txt -n {paste(n, collapse=',')} -o results") %>% system()
  res <- data.table::fread("results_post.txt")
  zs <- lapply(1:nid, function(i)
  {
    out <- list()
    out[[paste0("ZSCORE.P", i)]] <- regiondata[[i]]$x$beta.exposure / regiondata[[i]]$x$se.exposure
    as_tibble(out)
  }) %>% bind_cols()
  stopifnot(all(res$SNP_ID == regiondata[[1]]$x$SNP))
  names(res)[1] <- c("SNP")
  res <- bind_cols(res, zs) %>%
    mutate(CHR = regiondata[[1]]$x$chr.exposure, POS=regiondata[[1]]$x$position.exposure)
  setwd(thisdir)
  return(res)
}

#' Plot MsCAVIAR results
#'
#' @param res Output from run_MsCAVIAR
#'
#' @export
#' @return plot
plot_MsCAVIAR <- function(res)
{
  res %>% pivot_longer(cols=c(ZSCORE.P1, ZSCORE.P2, Prob_in_pCausalSet)) %>%
    ggplot(., aes(x=POS, y=value)) +
    geom_point(aes(colour=name)) +
    facet_grid(name ~ ., scale="free_y")
}

