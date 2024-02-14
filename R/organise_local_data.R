#' CAMERA_local class
#'
#' @description
#' A simple wrapper function for importing data from local files for use with the CAMERA class.
#' @export
CAMERA_local <- R6::R6Class("CAMERA_local", list(
    metadata = NULL,
    ld_ref = NULL,
    mc.cores = NULL,
    plink_bin = NULL,
    minmaf = NULL,
    pthresh = NULL,
    instrument_raw = NULL,
    instrument_outcome = NULL,
    instrument_regions = NULL,
    instrument_outcome_regions = NULL,

    # Methods
    #' @description
    #' Create a new dataset and initialise an R interface
    #' @param metadata Data frame with information about the data. One row per dataset. See details for info on columns
    #' @param ld_ref Data frame with two columns - pop = population (referencing the pop values in metadata), bfile = path to plink file for that reference
    #' @param plink_bin Location of executable plink (ver.1.90 is recommended)
    #' @param radius Genomic window size to extract SNPs
    #' @param clump_pop Reference population for clumping
    #' @param pthresh P-value threshold for instrument inclusion
    #' @param minmaf Minimum allelel frequency per dataset
    initialize = function(metadata, ld_ref, plink_bin, mc.cores=1, radius = 25000, pthresh = 5e-8, minmaf=0.01) {
        self$metadata <- metadata
        self$plink_bin <- plink_bin
        self$radius <- radius
        self$mc.cores <- mc.cores
        self$pthresh <- pthresh
        self$minmaf <- minmaf
    },

    standardise = function(d, ea_col="ea", oa_col="oa", beta_col="beta", eaf_col="eaf", chr_col="chr", pos_col="pos", vid_col="vid") {
        toflip <- d[[ea_col]] > d[[oa_col]]
        d[[eaf_col]][toflip] <- 1 - d[[eaf_col]][toflip]
        d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
        temp <- d[[oa_col]][toflip]
        d[[oa_col]][toflip] <- d[[ea_col]][toflip]
        d[[ea_col]][toflip] <- temp
        d[[vid_col]] <- paste0(d[[chr_col]], ":", d[[pos_col]], "_", d[[ea_col]], "_", d[[oa_col]])
        d
    },

    read_file = function(m, minmaf=0.01) {
        stopifnot(nrow(m) == 1)
        stopifnot(file.exists(m$fn))
        a <- data.table::fread(m$fn)
        message("Read ", nrow(a), " rows")
        b <- tibble(
            chr = a[[m$chr_col]],
            pos = as.numeric(a[[m$pos_col]]),
            eaf = as.numeric(a[[m$eaf_col]]),
            beta = as.numeric(a[[m$beta_col]]),
            se = as.numeric(a[[m$se_col]]),
            pval = as.numeric(a[[m$pval_col]]),
            ea = a[[m$ea_col]],
            oa = a[[m$oa_col]]
        ) %>% 
        filter(eaf > minmaf & eaf < (1-minmaf)) %>%
        standardise()
        return(b)
    },

    pool_tophits = function(rawdat, tophits, metadata, radius = 250000, pthresh = 5e-8, mc.cores = 10) {
        regions <- GRanges(
            seqnames = tophits$chr,
            ranges = IRanges(start=tophits$pos-radius, end=tophits$pos+radius),
            vid=tophits$vid, 
            pop=tophits$pop,
            trait=tophits$trait
        )
        region_list <- lapply(unique(tophits$trait), \(tr) {
            temp <- reduce(subset(regions, trait == tr))
            temp$trait <- tr
            temp
        })

        region_extract <- lapply(1:length(region_list), \(tr) {
            region <- region_list[[tr]]
            mclapply(1:length(region), \(i) {
                message(i, " of ", length(region))
                a <- lapply(1:nrow(metadata), \(j) {
                    subset(rawdat[[j]], chr == as.character(seqnames(region)[i]) & pos <= end(region)[i] & pos >= start(region)[i]) %>%
                        mutate(trait = metadata$trait[j], pop = metadata$pop[j], id = metadata$id[j])
                }) %>% bind_rows()
            })
        })

        pool <- lapply(1:length(region_extract), \(tr) {
            region <- region_extract[[tr]]
            mclapply(1:length(region), \(i) {
                target_trait <- region_list[[tr]]$trait[1]
                a <- region[[i]]
                k <- a %>% group_by(vid) %>% summarise(nstudies=n())
                a <- left_join(a, k, by="vid")
                k <- a %>% filter(trait == target_trait) %>%
                    group_by(nstudies) %>% 
                    summarise(minp = min(pval)) %>% 
                    filter(minp < pthresh)
                a <- subset(a, nstudies %in% k$nstudies)
                k <- subset(a, trait == target_trait) %>% 
                    mutate(z = abs(beta)/se) %>%
                    {subset(., z==max(z))$vid[1]}
                a <- subset(a, vid == k) %>% mutate(target_trait=target_trait)
                return(a)
            }, mc.cores=10) %>% 
                bind_rows()
        }) %>% bind_rows()

        region_list <- lapply(region_list, as_tibble)

        out <- list(region_list=region_list, region_extract=region_extract, tophit_pool=pool)
        return(out)
    },

    organise_data = function(metadata=self$metadata, plink_bin=self$plink_bin, ld_ref=self$ld_ref, pthresh=self$pthresh, minmaf = self$minmaf, radius = self$radius, mc.cores = self$mc.cores) {
        # read in data

        if(is.null(rawdat)) {
            rawdat <- mclapply(1:nrow(metadata), \(i) read_file(metadata[i,]))
        }

        # get top hits for each
        tophits <- lapply(1:nrow(metadata), \(i) {
            print(i)
            x <- rawdat[[i]] %>% 
                filter(pval < pthresh) %>%
                mutate(rsid = vid)
            if(nrow(x) > 1) {
                ieugwasr::ld_clump(x, plink_bin=plink_bin, bfile=subset(ld_ref, pop == metadata$pop[i])$bfile) %>%
                    select(-c(rsid)) %>%
                    mutate(pop=metadata$pop[i], trait=metadata$trait[i])
            } else {
                NULL
            }
        }) %>% bind_rows()

        # Get +-250kb region for every tophit
        # Get union of regions
        # Extract regions from each trait
        # Keep SNPs that have at least one GWAS sig and present in all
        # Clump to get tophits
        out <- pool_tophits(rawdat, tophits, metadata, radius = radius, pthresh = pthresh, mc.cores = mc.cores)
        return(out)
    },

    fixed_effects_meta_analysis_fast = function(beta_mat, se_mat) {
        w <- 1 / se_mat^2
        beta <- rowSums(beta_mat * w, na.rm=TRUE) / rowSums(w, na.rm=TRUE)
        se <- sqrt(1 / rowSums(w, na.rm=TRUE))
        pval <- pnorm(abs(beta / se), lower.tail = FALSE)
        return(pval)
    },

    organise = function() {
        metadata <- self$metadata
        ld_ref <- self$ld_ref
        exposure_trait <- subset(metadata, what = "exposure")$trait[1]
        outcome_trait <- subset(metadata, what = "outcome")$trait[1]
        
        # Read exposure in once
        metadata_exp <- subset(metadata, what == "exposure")
        rawdat <- mclapply(1:nrow(metadata_exp), \(i) read_file(metadata_exp[i,]), mc.cores=self$mc.cores)

        out <- subset(metadata, what == "outcome")$trait %>%
            unique() %>% {
            lapply(., \(x) {
                message(x)

                temp <- subset(metadata, trait == x)
                rawdat_this <- mclapply(1:nrow(temp), \(i) read_file(temp[i,]))
                rawdat_this <- c(rawdat, rawdat_this)

                a <- organise_data(subset(metadata, trait %in% c(exposure_trait, x)), self$plink_bin, self$ld_ref, self$mc.cores, rawdat=rawdat_this)
                return(a)
            })}

        o <- out[[1]]
        inst <- unique(subset(o$tophit_pool, target_trait == exposure_trait)$vid)
        inst_o <- unique(subset(o$tophit_pool, target_trait == outcome_trait)$vid)

        names(o$region_extract[[1]]) <- inst
        names(o$region_extract[[2]]) <- inst_o

        instrument_raw <- o$tophit_pool %>% filter(target_trait == exposure_trait & trait == exposure_trait) %>% rename(position="pos", nea="oa", p="pval", rsid="vid")
        instrument_outcome <- subset(o$tophit_pool, trait == outcome_trait & target_trait == exposure_trait & vid %in% instrument_raw$rsid) %>% rename(position="pos", nea="oa", p="pval", rsid="vid")

        # restrict regions to common snps

        instrument_raw
        instrument_regions <- lapply(unique(instrument_raw$rsid), \(x) {
            a <- o$region_extract[[1]][[x]] %>% 
                filter(trait == exposure_trait) %>% 
                rename(position="pos", nea="oa", p="pval", rsid="vid") %>%
                group_by(pop) %>% 
                group_split() %>% as.list()
            names(a) <- sapply(a, \(z) z$id[1])
            a
        })

        instrument_outcome_regions <- lapply(unique(instrument_raw$rsid), \(x) {
            a <- o$region_extract[[1]][[x]] %>% 
                filter(trait == outcome_trait) %>% 
                rename(position="pos", nea="oa", p="pval", rsid="vid") %>%
                group_by(pop) %>% 
                group_split() %>% as.list()
            names(a) <- sapply(a, \(z) z$id[1])
            a
        })

        names(instrument_regions) <- unique(instrument_raw$rsid)
        names(instrument_outcome_regions) <- unique(instrument_raw$rsid)

        self$instrument_regions <- instrument_regions
        self$instrument_outcome_regions <- instrument_outcome_regions
        self$instrument_raw <- instrument_raw
        self$instrument_outcome <- instrument_outcome
    }
))

