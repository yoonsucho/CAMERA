#' R6 class for CAMERA
#'
#' @description
#' A simple wrapper function.
#' Using a summary set, identify set of instruments for the traits, and peform SEM MR to test the association across the population.
#' @export

CAMERA <- R6::R6Class("CAMERA", list(
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
  mscaviar_results = NULL,
  expected_replications = NULL,
  instrument_region_zscores = NULL,
  instrument_maxz = NULL,
  instrument_susie = NULL,
  instrument_paintor= NULL,
  instrument_mscaviar = NULL,
  harmonised_data_check=NULL,
  standardised_instrument_raw=NULL,
  standardised_instrument_maxz=NULL,
  standardised_instrument_susie=NULL,
  standardised_instrument_paintor=NULL,
  standardised_instrument_mscaviar=NULL,
  standardised_outcome=NULL,
  instrument_specificity = NULL,
  instrument_specificity_summary = NULL,
  instrument_outcome = NULL,
  harmonised_dat_sem = NULL,
  sem_result = NULL,
  pleiotropic_snps = NULL,
  pleiotropy_dat = NULL,
  instrument_pleiotropy_summary = NULL,

# for convenience can migrate the results from a previous CAMERA into this one
#' @description
#' Migrate the results from a previous CAMERA
  import = function(x)
  {
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

#' @description
#'  Evaluate instrument associations between populations. The function checks whether 1) the instruments (chosen summary statistics IDs) can be used for the further steps and 2) the units are comparable between populations.
#' @param ids ID for the exposure or the outcome. Default is x$exposure_ids.
  check_phenotypes = function(ids=self$exposure_ids)
  {
    o <- lapply(ids, function(i) {tryCatch(
      {
        suppressMessages(exp <- unique(TwoSampleMR::extract_instruments(outcomes=i)))
        other_ids <- ids[!ids %in% i]

        o <- lapply(other_ids, function(j)
          {
            suppressMessages(out <- TwoSampleMR::extract_outcome_data(snps=exp$SNP, outcomes=j))
            suppressMessages(d <- TwoSampleMR::harmonise_data(exp, out))

            res <- suppressMessages(TwoSampleMR::mr(d, method="mr_ivw")) %>%
                     {tibble::tibble(Reference=i, Replication=j, nsnp=.$nsnp, agreement=.$b, se=.$se, pval=.$pval)}

            message(paste0("Instrument associations between ", i, " and ", j, " is ", round(res$agreement, 3), "; NSNP=", res$nsnp))

            suppressMessages(
              t <- d %>% data.frame() %>%
                  TwoSampleMR::add_metadata(., cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")) %>%
                  dplyr::summarise(unit_ref = .$units.exposure[1], unit_rep = .$units.outcome[1])
            )

            if(any(is.na(t))){
               message("Unit information is missing: See vignettes")
               print(t)
              }

            if(!any(is.na(t)) & t$unit_ref!=t$unit_rep){
                message("Units for the beta are different across the populations: Run x$standardise_data()")
                print(t)
              }
            return(res)
          })
        return(o %>% dplyr::bind_rows())
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})})
    print(o %>% dplyr::bind_rows())
  },

#' @description
#'  Identifies the instruments for the exposure
#' @param exposure_ids ID for the exposure. Default is x$exposure_ids.
  extract_instruments = function(exposure_ids=self$exposure_ids, ...)
  {
    # Use MVMR method to do the initial extraction
    # It gets the tophits in each exposure
    # Then randomly chooses one SNP per region to keep using clumping
    suppressMessages(instrument_raw <- TwoSampleMR::mv_extract_exposures(exposure_ids, ...))
    # Add chromosome and position
    suppressMessages(instrument_raw <- TwoSampleMR::add_metadata(instrument_raw, cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")))
    suppressMessages(instrument_raw <- ieugwasr::variants_rsid(unique(instrument_raw$SNP)) %>%
                         dplyr::select(SNP=query, chr, position=pos) %>%
                         dplyr::inner_join(., instrument_raw, by="SNP") %>%
                         dplyr::arrange(id.exposure, chr, position))
    # Arrange to be in order of exposure_ids
    # Rename columns
    instrument_raw <- lapply(self$exposure_ids, function(id) {
      subset(instrument_raw, id.exposure==id)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::select(rsid=SNP, chr, position, id=id.exposure, beta=beta.exposure, se=se.exposure, p=pval.exposure, ea=effect_allele.exposure, nea=other_allele.exposure, eaf=eaf.exposure, units=units.exposure, samplesize=contains("size")) %>%
      dplyr::mutate(method="raw") %>% as.data.frame
    # Check if top hits are significant in both populations
    t <- instrument_raw %>% dplyr::group_by(id) %>%
          dplyr::summarise(sum(p < 5e-8))
    id <- list()
    id <- t$id[t$`sum(p < 5e-08)` < 1]
    if(length(id) > 0){
      message(paste0("Caution: No SNPs reached genome-wide significance threshold for the trait in ", id))
    }
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
    self$instrument_regions <- lapply(regions, function(r) {tryCatch({
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
     } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
    names(self$instrument_regions) <- unique(instrument_raw$rsid)
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
    instrument_regions <- instrument_regions[lengths(instrument_regions) != 0]
    temp <- lapply(instrument_regions, function(r) {
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
          dplyr::mutate(original=r)
        instrument_regions[[r]] %>%
          lapply(., function(id){
            subset(id, rsid == o$rsid) %>%
              dplyr::mutate(original_rsid = r, zsum=o$zsum)
          }) %>% dplyr::bind_rows()
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(id, chr, position)
      t <- instrument_raw %>% dplyr::group_by(id) %>% dplyr::filter(dplyr::row_number()==1) %>%
           dplyr::select(id, units, samplesize)
      self$instrument_maxz  <- dplyr::left_join(self$instrument_maxz, t, by = "id") %>% as.data.frame()
    self$instrument_maxz <- lapply(self$exposure_ids, function(i)
    {
      subset(self$instrument_maxz, id == i) %>% dplyr::distinct(., rsid, .keep_all = TRUE)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(method = "maxz")
  },

#' @description
#' insert
#' @param region insert
#' @param instrument_region_zscores insert
#' @param instruments insert
#' @param comparison insert
  plot_regional_instruments = function(instrument_region_zscores=self$instrument_region_zscores, instruments=self$instrument_raw, region=1:min(10, nrow(instruments)), comparison=FALSE)
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

      if(comparison==FALSE){
            a$selected <- a$rsid %in% instruments$rsid
            ggplot2::ggplot(a, ggplot2::aes(x=position, y=value)) +
              ggplot2::geom_point(ggplot2::aes(colour=name)) +
              ggplot2::geom_point(data=subset(a, selected), , colour="black") +
              ggplot2::facet_grid(name ~ original_rsid, scale="free")
      }

      if(comparison==TRUE){
            a$selected_raw <- a$rsid %in% self$instrument_raw$rsid
            a$selected_maxz <- a$rsid %in% self$instrument_maxz$rsid
            a$selected_susie <- a$rsid %in% self$instrument_susie$rsid
            a$selected_paintor <- a$rsid %in% self$instrument_paintor$rsid
            ggplot2::ggplot(a, ggplot2::aes(x=position, y=value)) +
              ggplot2::geom_point(ggplot2::aes(colour=name)) +
              ggplot2::geom_point(data=subset(a, selected_raw), , colour="black") +
              ggplot2::geom_point(data=subset(a, selected_maxz), , colour="purple", shape=15) +
              ggplot2::geom_point(data=subset(a, selected_susie), , colour="orange", shape=17) +
              ggplot2::geom_point(data=subset(a, selected_paintor), , colour="brown", shape=18) +
              ggplot2::facet_grid(name ~ original_rsid, scale="free")
      }
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
    ld_matrices <- lapply(regions, function(r) {tryCatch(
    {
      d <- instrument_regions[[r]]
      exp <- self$exposure_ids
      o <- lapply(1:length(self$exposure_ids), function(i)
      {
        ld <- ieugwasr::ld_matrix(d[[exp[i]]]$rsid, pop=pops[i], bfile=bfiles[i], plink=self$plink, with_alleles=TRUE)
        code1 <- paste0(d[[exp[i]]]$rsid, "_", d[[exp[i]]]$ea, "_", d[[exp[i]]]$nea)
        code2 <- paste0(d[[exp[i]]]$rsid, "_", d[[exp[i]]]$nea, "_", d[[exp[i]]]$ea)
        rem_index <- ! (code1 %in% colnames(ld) | code2 %in% colnames(ld))
        flip_index <- ! colnames(ld) %in% code1
        if(any(rem_index))
        {
          rem <- d[[exp[i]]]$rsid[rem_index]
        }
        if(!any(rem_index))
        {
          rem <- NULL
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
      names(o) <- exp
      for(i in exp)
      {
        instrument_regions[[1]][[i]]$ld_unavailable <- instrument_regions[[1]][[i]]$rsid %in% rem
      }
      o <- lapply(o, function(ld)
      {
        rs <- strsplit(colnames(ld), "_") %>% sapply(., function(x) x[1])
        ld[! rs %in% rem, ! rs %in% rem]
      })
      return(o)
    } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
    names(ld_matrices) <- regions
    self$ld_matrices <- ld_matrices
    invisible(self)
  },

#' @description
#' insert
#' @param instrument insert
#' @param ld insert
  replication_evaluation = function(instrument=self$instrument_raw, ld=self$ld_matrices){
    instrument_pop <- instrument %>% dplyr::group_split(id)

    rep <- instrument_pop[[1]] %>% dplyr::select(rsid)
    rep$sign <- sign(instrument_pop[[1]]$beta)==sign(instrument_pop[[2]]$beta)
    rep$sig <- instrument_pop[[1]]$p<5e-8 & instrument_pop[[1]]$p<5e-8

    fam1 <- read.table(paste0(self$bfiles[[1]], ".fam"), header = FALSE)
    n1 <- length(unique(fam1[[1]]))

    ldsc_pop1 <- lapply(1:length(ld), function(i){
      ld <- ld[[i]][[1]]
      r2 <- ((n1-1) / (n1-2) * (ld^2)) - (1/(n1-2))
      l <- mean(r2[lower.tri(r2)]^2, diag=FALSE)
      return(l)
    }) %>% unlist()

    fam2 <- read.table(paste0(self$bfiles[[2]], ".fam"), header = FALSE)
    n2 <- length(unique(fam2[[1]]))

    ldsc_pop2 <- lapply(1:length(ld), function(i){
      ld <- ld[[i]][[2]]
      r2 <- ((n2-1) / (n2-2) * (ld^2)) - (1/(n2-2))
      l <- mean(r2[lower.tri(r2)]^2, diag=FALSE)
      return(l)
    }) %>% unlist()

    rep$delta_ld <- ldsc_pop1 - ldsc_pop2

    instrument_pop <- lapply(1:length(self$exposure_ids), function(i){
                               instrument_pop[[i]] %>% dplyr::mutate(maf = dplyr::if_else(.$eaf>0.5, (1-.$eaf), .$eaf, NA_real_))
                            })

    maf_pop1 <- instrument_pop[[1]]$maf * (1-instrument_pop[[1]]$maf)
    maf_pop2 <- instrument_pop[[2]]$maf * (1-instrument_pop[[2]]$maf)

    rep$delta_maf <- maf_pop1 - maf_pop2

    res <- list()
    res[[1]] <- summary(glm(sign ~ delta_ld, data = rep, family = binomial(link = "logit")))$coefficients
    res[[2]] <- summary(glm(sig ~ delta_maf, data = rep, family = binomial(link = "logit")))$coefficients
    res[[3]] <- summary(glm(sign ~ delta_ld + delta_maf, data = rep, family = binomial(link = "logit")))$coefficients
    res[[4]] <- summary(glm(sig ~ delta_ld + delta_maf, data = rep, family = binomial(link = "logit")))$coefficients
    names(res)[1] <- c("replicated_sign_model1")
    names(res)[2] <- c("replicated_sig_model1")
    names(res)[3] <- c("replicated_sign_model2")
    names(res)[4] <- c("replicated_sig_model2")
    return(list(res))
  },

  # for each instrument region + ld matrix we can perform susie finemapping
  # do this independently in each population
  # find credible sets that overlap - to try to determine the best SNP in a region to be used as instrument
#' @description
#' insert
#' @param dat insert
#' @param ld insert
 susie_finemap_regions = function(dat=self$instrument_regions, ld=self$ld_matrices)
  {
   resnps <- names(ld)
   exp <- self$exposure_ids
   susie <- lapply(1:length(resnps), function(r) {tryCatch({
        snp <- rownames(ld[[r]][[exp[1]]])[rownames(ld[[r]][[exp[1]]]) %in% rownames(ld[[r]][[exp[2]]])]
        snp <- strsplit(snp, "_") %>% sapply(., function(x) x[1])

        d <- lapply(1:length(exp), function(i){
                    index <- which(dat[[r]][[i]]$rsid  %in% snp)
                    d <- dat[[r]][[i]][index, ]
                    return(d)
              })

        su <- lapply(1:length(exp), function(i){
                    message("Running susie for ", exp[i])
                    susieR::susie_rss(z=d[[i]]$beta/d[[i]]$se, R=ld[[r]][[i]])
              })

        message("Finding overlaps [", r, "/", length(resnps), "]")
        suo <- private$susie_overlaps(su[[1]], su[[2]])

        out <- list(
                chr=d[[1]]$chr,
                position=d[[1]]$position,
                radius=self$radius,
                a1=d[[1]],
                a2=d[[2]],
                su1=su[[1]],
                su2=su[[2]],
                suo=suo
          )

     null <- c(is.null(su[[1]]$sets$cs), is.null(su[[2]]$sets$cs))
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
         temp <- dplyr::inner_join(d[[1]], d[[2]], by="rsid") %>%
           dplyr::mutate(pvalrank = rank(p.x) / nrow(d[[1]]) + rank(p.y) / nrow(d[[1]])) %>%
           dplyr::arrange(pvalrank)
         out$bestsnp <- temp$rsid[1]
       } else {
         out$cs_overlap <- TRUE
         out$bestsnp <- d[[1]]$rsid[suo$v[1]]
       }
     }

     if(out$type == "pop1")
     {
       out$cs_overlap <- FALSE
       out$bestsnp <- d[[1]] %>% dplyr::arrange(p) %>% {.$rsid[1]}
     }
     if(out$type == "pop2")
     {
       out$cs_overlap <- FALSE
       out$bestsnp <- d[[2]] %>% dplyr::arrange(p) %>% {.$rsid[1]}
     }
     return(out)
    } , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
   })

   names(susie) <- names(self$instrument_regions)
   self$susie_results <- susie

   #susie <- susie[!sapply(susie, is.null)]
   o <- unique(lapply(resnps, function(r) {susie[[r]]$bestsnp}) %>% unlist())
   instrument_susie <- lapply(resnps, function(r){tryCatch({
     lapply(exp, function(id){
          dat[[r]][[id]] %>% subset(., rsid %in% o) %>%
                            dplyr::bind_rows() %>%
                            dplyr::arrange(id, chr, position)})
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }) %>% dplyr::bind_rows() %>% as.data.frame()

   t <- self$instrument_raw %>%
          dplyr::group_by(id) %>% dplyr::filter(dplyr::row_number()==1) %>% dplyr::select(id, units, samplesize)

   instrument_susie <- dplyr::left_join(instrument_susie, t, by = "id") %>% as.data.frame()
   instrument_susie <- lapply(exp, function(i){subset(instrument_susie, id == i) %>% dplyr::distinct(., rsid, .keep_all = TRUE)}) %>% dplyr::bind_rows() %>%
                       dplyr::mutate(method = "susie")
   self$instrument_susie <- instrument_susie
   invisible(self)
  },

  # PAINTOR allows finemapping jointly across multiple populations
  # returns a posterior probability of inclusion for each SNP
  # Choose the variant with highest posterior probability and associations in each exposure
  #' @description insert
  #' @param region Output from extract_regional_data
  #' @param ld insert
  #' @param PAINTOR Path to PAINTOR executable. Default="PAINTOR"
  #' @param workdir Working directory. Default=tempdir()
  paintor_finemap_regions = function(region=self$instrument_regions, ld=self$ld_matrices, PAINTOR="PAINTOR", workdir=tempdir())
  {
    id <- self$exposure_ids
    nid <- length(region)
    snps <- lapply(1:nid, function(i) {region[[i]][[1]]$rsid})
    zs <- lapply(1:nid, function(i)
    {
      o <- list()
      lapply(1:length(id), function(id){
        o[[paste0("ZSCORE.P", id)]]<- region[[i]][[id]]$beta / region[[i]][[id]]$se
        return(tibble::as_tibble(o))
      }) %>% dplyr::bind_cols()
    })

    locus <- lapply(1:nid, function(i){ tryCatch({
      l <- list()
      l <- tibble::tibble(CHR = region[[i]][[1]]$chr, POS = region[[i]][[1]]$position, RSID = region[[i]][[1]]$rsid)
      l <- l %>% dplyr::bind_cols(., zs[[i]])
      ldsnp <- strsplit(rownames(ld[[i]][[1]]), "_") %>% sapply(., function(x) x[1])
      snp <- l$RSID %in% ldsnp
      index <- which(l$RSID %in% ldsnp)
      d <- l[index, ]
      return(d)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })

    anno <- lapply(1:nid, function(i){ tryCatch({
      tibble::tibble(null=rep(1, nrow(locus[[i]])))
       }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      })

    #write.files
    lapply(1:nid, function(i) {tryCatch({
      write.table(locus[i][[1]], file=file.path(workdir, paste0("Locus", i)), row=F, col=T, qu=F)
      write.table(anno[[i]], file=file.path(workdir, paste0("Locus", i, ".annotations")), row=F, col=T, qu=F)
      write.table(paste0("Locus", i), file=file.path(workdir, paste0("input.files", i)), row=F, col=F, qu=F)

      lapply(1:length(id), function(id)
      {
        write.table(ld[[i]][[id]], file=file.path(workdir, paste0("Locus", i, ".LD", id)), row=FALSE, col=FALSE, qu=FALSE)
      })
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })

    LDname <- paste0("LD", 1:length(id), collapse=',')
    Zhead <- paste(names(zs[[1]]), collapse=',')

    wd <- getwd()
    setwd(workdir)
    lapply(1:nid, function(i) { tryCatch({
      system(glue::glue("{PAINTOR} -input input.files{i} -Zhead {Zhead} -LDname {LDname} -in {workdir}/ -out {workdir}/ -mcmc -annotations null"))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })

    res <- lapply(1:nid, function(i) { tryCatch({
      data.table::fread(file.path(workdir, paste0("Locus", i, ".results"))) %>% dplyr::mutate(ZSCORE.SUM = ZSCORE.P1 + ZSCORE.P2) %>%
        dplyr::rename(., chr=CHR, position=POS, rsid=RSID)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
    names(res) <- names(self$instrument_regions)
    unlink(workdir)
    setwd(wd)
    self$paintor_results <- res

    #res <- res[!sapply(res, is.null)]
    bestsnp <- lapply(1:length(res), function(r){ tryCatch({
      if(sum(res[[r]]$Posterior_Prob)==0) {NULL}
      else {res[[r]] %>% dplyr::arrange(desc(Posterior_Prob)) %>% {.$rsid[1]}}
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })

    instrument_paintor <- lapply(1:nid, function(r){ tryCatch({
      lapply(1:length(id), function(id){
        region[[r]][[id]] %>% subset(., rsid %in% bestsnp) %>%
          dplyr::bind_rows() %>%
          dplyr::arrange(id, chr, position)
      })}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}) %>% dplyr::bind_rows() %>% as.data.frame()

    t <- self$instrument_raw %>% dplyr::group_by(id) %>% dplyr::filter(dplyr::row_number()==1) %>%
      dplyr::select(id, units, samplesize)
    instrument_paintor <- dplyr::left_join(instrument_paintor, t, by = "id") %>% as.data.frame()
    instrument_paintor <- lapply(id, function(i) {tryCatch({
      subset(instrument_paintor, id == i) %>% dplyr::distinct(., rsid, .keep_all = TRUE)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})}) %>% dplyr::bind_rows() %>%
      dplyr::mutate(method = "paintor")
    self$instrument_paintor <- instrument_paintor
    invisible(self)
  },


  #' @description Analyse regional data with MsCAVIAR
  #' @param region Output from extract_regional_data
  #' @param ld insert
  #' @param MsCAVIAR Path to MsCAVIAR executable. Default="MsCAVIAR"
  #' @param workdir Working directory. Default=tempdir()
  MsCAVIAR_finemap_regions = function(region=self$instrument_regions, ld=self$ld_matrices, MsCAVIAR="MsCAVIAR", workdir=tempdir())
  {
    id <- self$exposure_ids
    nid <- length(region)
    snps <- lapply(1:nid, function(i) {region[[i]][[1]]$rsid})
    n <- paste0(unique(self$instrument_raw$samplesize), collapse = ",")

    zs <- lapply(1:nid, function(i)
                {lapply(1:length(id), function(id)
                  {
                    zs <- region[[i]][[id]] %>%
                              dplyr::mutate(z=beta/se) %>%
                              dplyr::select(rsid, z)
                    names(zs)[2] <- c(paste0("ZSCORE.P", id))
                    ldsnp <- strsplit(rownames(ld[[i]][[1]]), "_") %>% sapply(., function(x) x[1])
                    snp <- zs$rsid %in% ldsnp
                    index <- which(zs$rsid %in% ldsnp)
                    zs <- zs[index, ]
                    write.table(zs, file=file.path(workdir, paste0("z_", i, "_", id, ".zscores")), row=F, col=F, qu=F)
                    return(tibble::as_tibble(zs))
                })
          })

    lapply(1:nid, function(i)
      {lapply(1:length(id), function(id)
        {
           write.table(ld[[i]][[id]], file=file.path(workdir, paste0("ld_", i, "_", id, ".ld")), row=FALSE, col=FALSE, qu=FALSE)
      })})

    thisdir <- getwd()
    setwd(workdir)

    lapply(1:nid, function(i)
      {
        writeLines(c(paste0("z_", i, "_", 1:length(id), ".zscores")), file(file.path(workdir, paste0("zfiles", i, ".txt"))))
        writeLines(c(paste0("ld_", i, "_", 1:length(id), ".ld")), file(file.path(workdir, paste0("ldfiles", i, ".txt"))))
      })

    lapply(1:nid, function(i)
              {
                system(glue::glue("{MsCAVIAR} -l ldfiles{i}.txt -z zfiles{i}.txt -n {n} -o results{i}"))
                message("Analysis [", i, "/", nid, "]")
              })

    res <- lapply(1:nid, function(i)
      {tryCatch({
        res <- data.table::fread(file.path(workdir, paste0("results", i, "_post.txt")))
        names(res)[1] <- c("rsid")
        names(res)[2] <- c("Posterior_Prob")
        return(res)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      })


    zs2 <- suppressMessages(lapply(1:nid, function(i)
                {
                  zs[[i]] %>% dplyr::bind_cols() %>%
                              dplyr::select(rsid=rsid...1, ZSCORE.P1, ZSCORE.P2)
                }))

    res <- suppressMessages(lapply(1:nid, function(i)
              {tryCatch({
                  res <- res[[i]] %>% dplyr::left_join(., zs2[[i]]) %>%
                              dplyr::left_join(., region[[i]][[1]]) %>%
                              dplyr::select(rsid, chr, position, ZSCORE.P1, ZSCORE.P2, Posterior_Prob)
                  return(res)
                  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }))
    names(res) <- names(self$instrument_regions)
    unlink(workdir)
    setwd(thisdir)

    self$mscaviar_results <- res

    bestsnp <- lapply(1:length(res), function(r)
      {
        if(sum(res[[r]]$Posterior_Prob)==0) {NULL}
        else {res[[r]] %>% dplyr::arrange(desc(Posterior_Prob)) %>% {.$rsid[1]}}
      })

    instrument_mscaviar <- lapply(1:nid, function(r){
        lapply(1:length(id), function(id){
          region[[r]][[id]] %>% subset(., rsid %in% bestsnp) %>%
            dplyr::bind_rows() %>%
            dplyr::arrange(id, chr, position)
        })}) %>% dplyr::bind_rows() %>% as.data.frame()

      t <- self$instrument_raw %>% dplyr::group_by(id) %>% dplyr::filter(dplyr::row_number()==1) %>%
        dplyr::select(id, units, samplesize)
      instrument_mscaviar <- dplyr::left_join(instrument_mscaviar, t, by = "id") %>% as.data.frame()
      instrument_mscaviar <- lapply(id, function(i) {subset(instrument_mscaviar, id == i)}) %>% dplyr::bind_rows() %>%
                             dplyr::relocate(rsid, chr, position, id, beta, se, p, ea, nea, units, samplesize) %>%
                             dplyr::mutate(method = "mscaviar")
      self$instrument_mscaviar <- instrument_mscaviar
      invisible(self)
  },

#' @description
#' Evaluate heterogeneity of identified instrument-trait associations across populations, where the instruments were identified using "Raw", "MaxZ", or fine-mapping (Susie, PAINTOR) methods.
#' @param instrument Intsruments for the exposure that are selected the provided methods (x$instrument_raw, x$instrument_maxz, x$instrument_susie, x$instrument_paintor). Default is x$instrument_raw.
#' @param alpha insert
#' @param method insert
 instrument_heterogeneity = function(instrument=self$instrument_raw, alpha="bonferroni", method="ivw")
  {
    if(alpha=="bonferroni")
    {
      alpha <- 0.05/nrow(instrument)
    }

    if(method=="ivw")
    {
     if(!any(names(instrument) %in% c("beta.outcome")))
     {
      o <- lapply(self$exposure_ids, function(i)
      {
        m <- subset(instrument, id == i & p < alpha)
        other_ids <- self$exposure_ids[!self$exposure_ids %in% i]

        o <- lapply(other_ids, function(j)
            {
              n <- subset(instrument, id == j & rsid %in% m$rsid)
              dat <-dplyr:: inner_join(m, n, by="rsid") %>%
                    dplyr::select(SNP=rsid, x=beta.x, y=beta.y, xse=se.x, yse=se.y, xp=p.x, yp=p.y)
              res <- suppressMessages(TwoSampleMR::mr_ivw(dat$x, dat$y, dat$xse, dat$yse)) %>%
                       {tibble::tibble(Reference=i, Replication=j, nsnp=length(unique(dat$SNP)), agreement=.$b, se=.$se, pval=.$pval, Q=.$Q, Q_pval=.$Q_pval, I2=((.$Q - length(dat$SNP))/.$Q))}
              return(res)
        })})
     }

     if(any(names(instrument) %in% c("beta.outcome")))
     {
      o <- lapply(self$outcome_ids, function(i)
      {
        m <- subset(instrument, id.outcome == i & pval.outcome < alpha)
        other_ids <- self$outcome_ids[!self$outcome_ids %in% i]

        o <- lapply(other_ids, function(j)
            {
              n <- subset(instrument, id.outcome == j & SNP %in% m$SNP)
              dat <-dplyr:: inner_join(m, n, by="SNP") %>%
                    dplyr::select(SNP=SNP, x=beta.outcome.x, y=beta.outcome.y, xse=se.outcome.x, yse=se.outcome.y, xp=pval.outcome.x, yp=pval.outcome.y)
              res <- suppressMessages(TwoSampleMR::mr_ivw(dat$x, dat$y, dat$xse, dat$yse)) %>%
                      {tibble::tibble(Reference=i, Replication=j, nsnp=length(unique(dat$SNP)), agreement=.$b, se=.$se, pval=.$pval, Q=.$Q, Q_pval=.$Q_pval, I2=((.$Q - length(dat$SNP))/.$Q))}
              return(res)
        })})
     }
    }

   if(method=="simple_mode")
   {
     if(!any(names(instrument) %in% c("beta.outcome")))
     {
       o <- lapply(self$exposure_ids, function(i)
       {
         m <- subset(instrument, id == i & p < alpha)
         other_ids <- self$exposure_ids[!self$exposure_ids %in% i]

         o <- lapply(other_ids, function(j)
         {
           n <- subset(instrument, id == j & rsid %in% m$rsid)
           dat <-dplyr:: inner_join(m, n, by="rsid") %>%
             dplyr::select(SNP=rsid, x=beta.x, y=beta.y, xse=se.x, yse=se.y, xp=p.x, yp=p.y)
           res <- suppressMessages(TwoSampleMR::mr_simple_mode_nome(dat$x, dat$y, dat$xse, dat$yse)) %>%
             {tibble::tibble(Reference=i, Replication=j, nsnp=length(unique(dat$SNP)), agreement=.$b, se=.$se, pval=.$pval)}
           return(res)
         })})
     }

     if(any(names(instrument) %in% c("beta.outcome")))
     {
       o <- lapply(self$outcome_ids, function(i)
       {
         m <- subset(instrument, id.outcome == i & pval.outcome < alpha)
         other_ids <- self$outcome_ids[!self$outcome_ids %in% i]

         o <- lapply(other_ids, function(j)
         {
           n <- subset(instrument, id.outcome == j & SNP %in% m$SNP)
           dat <-dplyr:: inner_join(m, n, by="SNP") %>%
             dplyr::select(SNP=SNP, x=beta.outcome.x, y=beta.outcome.y, xse=se.outcome.x, yse=se.outcome.y, xp=pval.outcome.x, yp=pval.outcome.y)
           res <- suppressMessages(TwoSampleMR::mr_simple_mode_nome(dat$x, dat$y, dat$xse, dat$yse)) %>%
             {tibble::tibble(Reference=i, Replication=j, nsnp=length(unique(dat$SNP)), agreement=.$b, se=.$se, pval=.$pval)}
           return(res)
         })})
     }

   }
    print(o %>% dplyr::bind_rows())
  },

#' @description
#' insert
#' @param dat insert
#' @param standardise_unit insert
#' @param standardise_scale insert
  standardise_data = function(dat=self$instrument_raw, standardise_unit=FALSE, standardise_scale=FALSE)
  {
    if(standardise_unit==TRUE)
    {
      if(!any(names(dat) %in% c("beta.outcome")))
      {
       exp <- dat
       d <- exp %>%
              dplyr::group_by(id) %>% dplyr::summarise(units = units[1])
     
       if(any(is.na(exp$eaf))){exp <- private$allele_frequency(dat = exp)}
     
       if(!any(d$units %in% c("log odds"))) {
           exp <- exp %>%
                       dplyr::group_by(id) %>%
                       dplyr::mutate(units = dplyr::na_if(units, "NA")) %>%
                       dplyr::mutate(units = replace(units, is.na(units), "temp")) %>%
                       dplyr::mutate(estimated_sd = mean(TwoSampleMR::estimate_trait_sd(beta, se, samplesize, eaf), na.rm=TRUE)) %>%
                       dplyr::mutate(estimated_sd = replace(estimated_sd, units=="SD", 1))}

       if(!any(is.na(exp$estimated_sd))) {
        stopifnot(!any(is.na(exp$estimated_sd)))
        exp$beta <- exp$beta / exp$estimated_sd
        exp$se <- exp$se / exp$estimated_sd
        exp$units <- "SD"

        exp <- exp %>% as.data.frame()

        if(any(exp$method[[1]] %in% c("raw"))){self$standardised_instrument_raw <- exp}
        if(any(exp$method[[1]] %in% c("maxz"))){self$standardised_instrument_maxz <- exp}
        if(any(exp$method[[1]] %in% c("susie"))){self$standardised_instrument_susie <- exp}
        if(any(exp$method[[1]] %in% c("paintor"))){self$standardised_instrument_paintor <- exp}
        if(any(exp$method[[1]] %in% c("mscaviar"))){self$standardised_instrument_mscaviar <- exp}}
     }

     if(any(names(dat) %in% c("beta.outcome")))
      {
        out <- dat
        d <- out %>%
             dplyr::group_by(id.outcome) %>% dplyr::summarise(units = units.outcome[1])

        if(any(is.na(out$eaf.outcome))){out <- private$allele_frequency(dat = out)}      

        if(!any(d$units %in% c("log odds"))){
          out <- out %>%
              dplyr::group_by(id.outcome) %>%
              dplyr::mutate(units.outcome = dplyr::na_if(units.outcome, "NA")) %>%
              dplyr::mutate(units.outcome = replace(units.outcome, is.na(units.outcome), "temp")) %>%
              dplyr::mutate(estimated_sd = mean(TwoSampleMR::estimate_trait_sd(beta.outcome, se.outcome, samplesize.outcome, eaf.outcome), na.rm=TRUE)) %>%
              dplyr::mutate(estimated_sd = replace(estimated_sd, units.outcome=="SD", 1))}

        if(any(is.na(out$estimated_sd))){ 
          stopifnot(!any(is.na(out$estimated_sd)))
          out$beta.outcome <- out$beta.outcome / out$estimated_sd
          out$se.outcome <- out$se.outcome / out$estimated_sd
          out$units.outcome <- "SD"
        }
        
        self$standardised_outcome <- out %>% as.data.frame()}
    }

    if(standardise_scale==TRUE)
    {
     if(!any(names(dat) %in% c("beta.outcome")))
      {
        stopifnot(!is.null(self$instrument_maxz))
        invisible(capture.output(scale <- self$instrument_heterogeneity(instrument=self$instrument_maxz, method="simple_mode")))
        bxx <- scale$agreement[1]
        if(!is.null(self[[paste0("standardised_instrument_", dat$method[[1]])]]))
        {
          dat <- self[[paste0("standardised_instrument_", dat$method[[1]])]]
        }
        exp <- dat %>%
            dplyr::mutate(original_beta = beta) %>%
            dplyr::mutate(original_se = se) %>%
            dplyr::mutate(beta = dplyr::case_when(id == self$exposure_ids[1] ~ beta * bxx, TRUE ~ beta)) %>%
            dplyr::mutate(se = dplyr::case_when(id == self$exposure_ids[1] ~ se * bxx, TRUE ~ se)) %>%
            as.data.frame()

        if(any(exp$method[[1]] %in% c("raw"))){self$standardised_instrument_raw <- exp}
        if(any(exp$method[[1]] %in% c("maxz"))){self$standardised_instrument_maxz <- exp}
        if(any(exp$method[[1]] %in% c("susie"))){self$standardised_instrument_susie <- exp}
        if(any(exp$method[[1]] %in% c("paintor"))){self$standardised_instrument_paintor <- exp}
        if(any(exp$method[[1]] %in% c("mscaviar"))){self$standardised_instrument_mscaviar <- exp}
      }

    if(any(names(dat) %in% c("beta.outcome")))
      {
        stopifnot(!is.null(self$instrument_outcome))

        if(!is.null(self$standardised_outcome))
        {
          dat <- self$standardised_outcome
        }

        oexp <- self$exposure_ids
        oraw <- self$instrument_raw
        omaxz <- self$instrument_maxz
        ore <- self$instrument_regions
        oz <- self$instrument_region_zscores

        self$exposure_ids <- self$outcome_ids
        suppressMessages(self$extract_instruments(exposure_ids=self$exposure_ids))
        suppressMessages(self$extract_instrument_regions(instrument_raw=self$instrument_raw, exposure_ids=self$exposure_ids))
        suppressMessages(self$scan_regional_instruments(instrument_raw=self$instrument_raw, instrument_regions=self$instrument_regions))
        invisible(capture.output(scale <- self$instrument_heterogeneity(instrument=self$instrument_maxz, method="simple_mode")))

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
#' @param winnerscurse insert
  estimate_instrument_specificity = function(instrument, alpha="bonferroni", winnerscurse=FALSE)
  {
    if(alpha=="bonferroni")
    {
      alpha <- 0.05/nrow(instrument)
    }
    o <- lapply(self$exposure_ids, function(i)
    {
      m <- subset(instrument, id == i & p < alpha)
      other_ids <- self$exposure_ids[!self$exposure_ids %in% i]

      if(winnerscurse==TRUE){
          wcm <- m %>% dplyr::select(rsid, beta, se) %>% dplyr::mutate(rsid = as.numeric(gsub("rs","", .$rsid)))
          wcl <- winnerscurse::cl_interval(summary_data=wcm, alpha = alpha, conf_level=0.95) %>%
                                  dplyr::mutate(rsid = sub("^", "rs",  .$rsid)) %>%
                                  dplyr::mutate(se.cl3 = (.$upper - .$lower) / 3.92) %>% dplyr::arrange(rsid)
          o <- lapply(other_ids, function(j)
          {
            n <- subset(instrument, id == j & rsid %in% wcl$rsid) %>% dplyr::arrange(rsid)
            o <- n %>%
              {private$prop_overlap(wcl$beta.cl3, .$beta, wcl$se.cl3, .$se, alpha)}
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
        })}

      if(winnerscurse==FALSE){
         o <- lapply(other_ids, function(j)
          {
            n <- subset(instrument, id == j & rsid %in% m$rsid)
            o <- n %>%
              {private$prop_overlap(m$beta, .$beta, m$se, .$se, alpha)}
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
      })}
      overall <- lapply(o, function(x) { x$res }) %>% dplyr::bind_rows()
      pervariant <- lapply(o, function(x) { x$variants }) %>% dplyr::bind_rows()
      return(list(overall=overall, pervariant=pervariant))
    })
    self$instrument_specificity_summary <- lapply(o, function(x) x$overall) %>% dplyr::bind_rows()
    self$instrument_specificity <- lapply(o, function(x) x$pervariant) %>% dplyr::bind_rows()
    return(self$instrument_specificity_summary)
  },

  # Once a set of instruments is defined then extract the outcome data
  # Could use
  # - raw results from extract_instruments
  # - maximised associations from extract_instrument_regions
  # - finemapped hits from susie_finemap_regions
  # - finemapped hits from mscaviar_finemap_regions
  #' @description
  #' insert
  #' @param exp insert
  #' @param p_exp insert
  make_outcome_data = function(exp=self$instrument_raw, p_exp=0.05/nrow(exp)){
      dx <- dplyr::inner_join(
                              subset(exp, id == self$exposure_ids[[1]]),
                              subset(exp, id == self$exposure_ids[[2]]),
                              by="rsid"
                              ) %>% dplyr::filter(p.x < p_exp)
      suppressMessages(out <- TwoSampleMR::extract_outcome_data(snps=dx$rsid, outcomes=self$outcome_ids))
      suppressMessages(out <- TwoSampleMR::add_metadata(out, cols = c("sample_size", "ncase", "ncontrol", "unit", "sd")) %>%
                              #dplyr::select(rsid=SNP, chr, position=pos, id=id.outcome, beta=beta.outcome, se=se.outcome, p=pval.exposure, ea=effect_allele.outcome, nea=other_allele.outcome, eaf=eaf.outcome, units=units.outcome, samplesize=contains("size")) %>%
                              as.data.frame()
                      )
      self$instrument_outcome <- out %>% dplyr::arrange(., chr, SNP)
      invisible(self)
  },

#' @description
#' insert
#' @param exp insert
#' @param out insert
  harmonised_dat = function(exp=self$instrument_raw, out=self$instrument_outcome){
      dx <- dplyr::inner_join(
                              subset(exp, id == self$exposure_ids[[1]]),
                              subset(exp, id == self$exposure_ids[[2]]),
                              by="rsid") %>%
            dplyr::select(SNP=rsid, x1=beta.x, x2=beta.y, xse1=se.x, xse2=se.y, p1=p.x, p2=p.y)
      dy <- dplyr::inner_join(
                              subset(out, id.outcome == self$outcome_ids[[1]]),
                              subset(out, id.outcome == self$outcome_ids[[2]]),
            by="SNP") %>%
          dplyr::select(SNP=SNP, y1=beta.outcome.x, y2=beta.outcome.y, yse1=se.outcome.x, yse2=se.outcome.y)
      dat <- dplyr::inner_join(dx, dy, by="SNP")
      self$harmonised_dat_sem <- dat
  },

# Generate harmonised dataset
# Perform basic SEM analysis of the joint estimates in multiple ancestries
#' @description
#' insert
#' @param harmonised_dat insert
  perform_basic_sem = function(harmonised_dat = self$harmonised_dat_sem) {
      d <- harmonised_dat %>%
           dplyr::mutate(r1 = y1/x1) %>%
           dplyr::mutate(r2 = y2/x2) %>%
           dplyr::mutate(w1 = sqrt(x1^2 / yse1^2)) %>%
           dplyr::mutate(w2 = sqrt(x2^2 / yse2^2)) %>%
           dplyr::mutate(o1 = r1 * w1) %>%
           dplyr::mutate(o2 = r2 * w2)

      out <- list()
      out$ivw1 <- TwoSampleMR::mr_ivw(d$x1, d$y1, d$xse1, d$yse1) %>%
                    {tibble::tibble(Methods="IVW", pop="1", nsnp=nrow(d), bivhat=.$b, se=.$se, pval=.$pval)}
      out$ivw2 <- TwoSampleMR::mr_ivw(d$x2, d$y2, d$xse2, d$yse2) %>%
                    {tibble::tibble(Methods="IVW", pop="2", nsnp=nrow(d), bivhat=.$b, se=.$se, pval=.$pval)}
      out$rm1 <- summary(lm(o1 ~ -1 + w1, data=d)) %>%
                    {tibble::tibble(Methods="RadialIVW", pop="1", nsnp=nrow(d), bivhat=.$coef[1,1], se=.$coef[1,2], pval=.$coef[1,4])}
      out$rm2 <- summary(lm(o2 ~ -1 + w2, data=d)) %>%
                    {tibble::tibble(Methods="RadialIVW", pop="2", nsnp=nrow(d), bivhat=.$coef[1,1], se=.$coef[1,2], pval=.$coef[1,4])}

      out$semA <- private$runsem('
                                 y1 ~ biv*x1
                                 y2 ~ biv*x2
                                 ', d, "UnweightedSEMa")[1, ] %>%
                                dplyr::mutate(pop=replace(pop, pop==1, "1=2"))

      out$semB <- private$runsem('
                                   y1 ~ biv_1*x1
                                   y2 ~ biv_2*x2
                                   ', d, "UnweightedSEMb")

      out$modC <- private$runsem('
                                 o1 ~ biv*w1
                                 o2 ~ biv*w2
                                 ', d, "RadialSEMa")[1, ] %>%
                                dplyr::mutate(pop=replace(pop, pop==1, "1=2"))

      out$modD <- private$runsem('
                                   o1 ~ biv_1*w1
                                   o2 ~ biv_2*w2
                                   ', d, "RadialSEMb")

      invisible(self$sem_result <- out %>% dplyr::bind_rows())
      print(self$sem_result)
  },


#' @description Return a list of outlying SNPs in each population 
#' @param harmonised_dat harmonised dataset generated using \code{x$harmonised_dat()}
  pleiotropy = function(harmonised_dat=self$harmonised_dat_sem){
    stopifnot(!is.null(harmonised_dat))

    sig <- harmonised_dat[harmonised_dat$p1<5*10^-8,]

    d1 <- RadialMR::format_radial(BXG=sig$x1, BYG=sig$y1, seBXG=sig$xse1, seBYG=sig$yse1, RSID=sig$SNP)
    invisible(capture.output(o1 <- RadialMR::ivw_radial(d1, alpha=1, weights=3)))
    o1$outliers$residuals <- summary(lm(BetaWj ~ -1 + Wj, o1$data))$residuals

    d2 <- RadialMR::format_radial(BXG=sig$x2, BYG=sig$y2, seBXG=sig$xse2, seBYG=sig$yse2, RSID=sig$SNP)
    invisible(capture.output(o2 <- RadialMR::ivw_radial(d2, alpha=1, weights=3)))
    o2$outliers$residuals <- summary(lm(BetaWj ~ -1 + Wj, o2$data))$residuals

    dat <- merge(o1$outliers, o2$outliers, by="SNP") %>%
           dplyr::mutate(
                          sigx=p.adjust(p.value.x, "fdr") < 0.05,
                          sigy=p.adjust(p.value.y, "fdr") < 0.05,
                          outlier = dplyr::case_when(
                            sigx & sigy ~ "Both",
                            sigx & ! sigy ~ "pop1",
                            ! sigx & sigy ~ "pop2",
                            ! sigx & ! sigy ~ "None")) %>%
           dplyr::mutate(outp1 = dplyr::if_else(outlier=="pop1", 1, dplyr::if_else(outlier=="both", 1, 0))) %>%
           dplyr::mutate(outp2 = dplyr::if_else(outlier=="pop2", 1, dplyr::if_else(outlier=="both", 1, 0))) 

    invisible(self$pleiotropy_dat <- dat %>% dplyr::select("SNP", "outlier", "outp1", "outp2"))

    out <- list()
    out$both <- subset(dat$SNP, dat$outlier=="Both")
    out$pop1 <- subset(dat$SNP, dat$outlier=="pop1")
    out$pop2 <- subset(dat$SNP, dat$outlier=="pop2")
    out$none <- subset(dat$SNP, dat$outlier=="None")

    invisible(self$pleiotropic_snps <- out)

    dat %>%
      ggplot2::ggplot(., ggplot2::aes(x=Q_statistic.x, y=Q_statistic.y)) +
      ggplot2::geom_point(ggplot2::aes(colour=outlier)) +
      #ggplot2::geom_smooth(method = "lm", se=FALSE, color="gray", alpha = .2, size = 0.2, formula = y ~ x) +
      ggplot2::scale_colour_brewer(type="qual") +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::ylab('Q statistics in pop2') +
      ggplot2::xlab('Q statistics in pop1')
  },


#' @description Estimate population specificity of pleiotric SNPs
#' @param harmonised_dat harmonised dataset generated using \code{x$harmonised_dat()}
#' @param sem_result MR-SEM result obtained by \code{x$perform_basic_sem()}
 pleiotropy_specificity = function(harmonised_dat=self$harmonised_dat_sem, sem_result=self$sem_result, pleioropy=self$pleiotropy_dat){
  stopifnot(!is.null(harmonised_dat))
  stopifnot(!is.null(sem_result))

  sig <- harmonised_dat[harmonised_dat$p1<5*10^-8,]

  if(sem_result$aic[5] - sem_result$aic[6] <=-2){
      d <- sig %>%
           dplyr::mutate(wald1=y1/x1, wald.se1= yse1/abs(x1), wald2=y2/x2, wald.se2= yse2/abs(x2),
                         ivw1=sem_result$bivhat[5], ivw.se1=sem_result$se[5], ivw2=sem_result$bivhat[5], ivw.se2=sem_result$se[5]) 
    }

  if(is.na(sem_result$aic[6])){
    message("Caution: SE is not properly estimated for model 2. The estimates from model 1 are used.")
      d <- sig %>%
           dplyr::mutate(wald1=y1/x1, wald.se1= yse1/abs(x1), wald2=y2/x2, wald.se2= yse2/abs(x2),
                         ivw1=sem_result$bivhat[5], ivw.se1=sem_result$se[5], ivw2=sem_result$bivhat[5], ivw.se2=sem_result$se[5]) 
    }               

  if(sem_result$aic[5] - sem_result$aic[6] > -2){
      d <- sig %>%
            dplyr::mutate(wald1=y1/x1, wald.se1= yse1/abs(x1), wald2=y2/x2, wald.se2= yse2/abs(x2),
                         ivw1=sem_result$bivhat[6], ivw.se1=sem_result$se[6], ivw2=sem_result$bivhat[7], ivw.se2=sem_result$se[7])
    } 

  pop1 <- list()
  pop2 <- list()
  for(i in 1:nrow(d))
      {
       pop1[[i]] <- private$bootstrap(d$wald1[i], d$wald.se1[i], d$ivw1[i], d$ivw.se1[i], nboot=1000)
       pop2[[i]] <- private$bootstrap(d$wald2[i], d$wald.se2[i], d$ivw2[i], d$ivw.se2[i], nboot=1000)
      }
  pop1 <- do.call("rbind", pop1) %>% as.data.frame %>% dplyr::select(pleio.p1=pleio, sd.p1=sd) 
  pop2 <- do.call("rbind", pop2) %>% as.data.frame %>% dplyr::select(pleio.p2=pleio, sd.p2=sd)
  
  d <- cbind(d, pop1, pop2) %>%
            merge(., pleioropy, by = "SNP")
  
  p1 <- d %>% dplyr::select("SNP", x="x1", xse="xse1", p="p1", y="y1", yse="yse1", wald="wald1", wald.se1="wald.se1", pleio="pleio.p1", sd="sd.p1", out="outp1") %>%
                  dplyr::mutate(id=self$exposure_ids[1])
  p2 <- d %>% dplyr::select("SNP", x="x2", xse="xse2", p="p2", y="y2", yse="yse2", wald="wald2", wald.se1="wald.se2", pleio="pleio.p2", sd="sd.p2", out="outp2") %>%
                  dplyr::mutate(id=self$exposure_ids[2]) 

  mer <- rbind(p1, p2)

  o <- lapply(self$exposure_ids, function(i)
       {
        m <- subset(mer, id == i & out==1)
        other_ids <- self$exposure_ids[!self$exposure_ids %in% i]
        o <- lapply(other_ids, function(j)
             {
              n <- subset(mer, id == j & SNP %in% m$SNP) %>% dplyr::arrange(SNP)
              o <- n %>%
                    {private$prop_overlap(m$pleio, .$pleio, sqrt(m$sd), sqrt(.$sd), alpha=0.05)}
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
              return(o)})
        overall <- lapply(o, function(x) { x$res }) %>% dplyr::bind_rows()
        pervariant <- lapply(o, function(x) { x$variants }) %>% dplyr::bind_rows()
        return(list(overall=overall, pervariant=pervariant))
        })

      self$instrument_pleiotropy_summary <- lapply(o, function(x) x$overall) %>% dplyr::bind_rows()
      print(self$instrument_pleiotropy_summary)
  }

))



