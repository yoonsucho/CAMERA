#' Fixed effects meta analysis vectorised across multiple SNPs
#'
#' Assumes effects across studies are all on the same scale
#'
#' @param beta_mat Matrix of betas - rows are SNPs, columns are studies
#' @param se_mat Matrix of SEs - rows are SNPs, columns are studies
#'
#' @return list of meta analysis betas and SEs
#' @export
fixed_effects_meta_analysis_fast <- function(beta_mat, se_mat) {
  w <- 1 / se_mat^2
  beta <- rowSums(beta_mat * w) / rowSums(w, na.rm=TRUE)
  se <- sqrt(1 / rowSums(w, na.rm=TRUE))
  z <- abs(beta / se)
  p <- pnorm(z, lower.tail = FALSE)
  nstudy <- apply(beta_mat, 1, \(x) sum(!is.na(x)))
  return(tibble(nstudy, p, z=z))
}

#' P-value based meta analysis
#'
#' Uses weighted Z scores following advice from https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x
#' Suggested weights are 1/se^2
#' However se is scale dependent and it would be ideal to avoid scale issues at this stage
#' So using instead calculate expected se based on n and af.
#' Assumes continuous traits in the way it uses n (i.e. not case control aware at the moment)
#'
#' @param beta_mat Matrix of betas - rows are SNPs, columns are studies
#' @param se_mat Matrix of SEs - rows are SNPs, columns are studies
#' @param n Vector of sample sizes for each
#' @param eaf_mat Matrix of allele frequencies - rows are SNPs, columns are studies
z_meta_analysis <- function(beta_mat, se_mat, n, eaf_mat) {
  z_mat <- abs(beta_mat) / se_mat
  w_mat <- t(t(sqrt(eaf_mat * (1 - eaf_mat) * 2)) * sqrt(n))
  zw <- rowSums(z_mat * w_mat, na.rm=TRUE) / sqrt(rowSums(w_mat^2, na.rm=TRUE))
  p <- pnorm(abs(zw), lower.tail = FALSE)
  nstudy <- apply(beta_mat, 1, \(x) sum(!is.na(x)))
  return(tibble(nstudy, p, z=abs(zw)))
}

#' Identify best variant for each region
#' 
#' Use the 
#' 
#' @param dat Output from `pleiotropy` - `pleiotropy_outliers`
#' 
#' @return plot
CAMERA$set("public", "fema_regional_instruments", function(method = "fema", instrument_regions = self$instrument_regions, instrument_raw = self$instrument_raw, n=self$exposure_metadata$sample_size) {
  
  stopifnot(method %in% c("fema", "zma"))
  if(method == "zma") {
    stopifnot(!is.null(n))
    stopifnot(all(!is.na(n)))
  }
  
  # remove duplicated ids from instrument_regions
  instrument_regions <- lapply(instrument_regions, \(x) {
    lapply(x, \(y) {
      subset(y, !duplicated(rsid))
    })
  })

  # remove regions that have no data
  rem <- lengths(instrument_regions) == 0
  if (any(rem)) {
    message(sum(rem), " regions have no data")
    instrument_regions <- instrument_regions[!rem]
  }

  # for every region, create mats, do scan
  d <- lapply(instrument_regions, \(x) {
    if(is.null(x)) return(NULL)
    if(nrow(x[[1]]) == 0) return(NULL)
      
      rsidintersect <- Reduce(intersect, lapply(x, \(r) r$rsid))
      x <- lapply(x, \(r) {
        r %>% dplyr::filter(rsid %in% rsidintersect) %>% dplyr::filter(!duplicated(rsid)) %>% dplyr::arrange(rsid)
      })
      
      d <- dplyr::select(x[[1]], chr, position, rsid, ea, nea, rsido, trait)

      if(method == "fema") {
        d1 <- fixed_effects_meta_analysis_fast(
          sapply(x, \(y) y$beta), 
          sapply(x, \(y) y$se)
        )
      } else {
        d1 <- z_meta_analysis(
          sapply(x, \(y) y$beta),
          sapply(x, \(y) y$se),
          n,
          sapply(x, \(y) y$eaf)
        )
      }
      d <- dplyr::bind_cols(d, d1)
  })

  # Keep best SNP from each region
  names(d) <- names(instrument_regions)
  d[sapply(d, is.null)] <- NULL
  d[sapply(d, \(d1) {nrow(d1) == 0})] <- NULL
  dsel <- lapply(d, \(x) {
    subset(x, z == max(x$z, na.rm=TRUE))[1,]
  }) %>% dplyr::bind_rows()
  dsel$region <- names(d)
  # Extract best SNPs from regions for each pop
  inst <- lapply(1:nrow(dsel), \(i) {
    lapply(instrument_regions[[dsel$region[i]]], \(p) {
      
      subset(p, rsid == dsel$rsid[i])
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(id=names(instrument_regions[[dsel$region[i]]]))
  }) %>% dplyr::bind_rows() %>%
    dplyr::filter(!duplicated(paste(id, rsid)))
  self$instrument_fema <- inst
  self$instrument_fema_regions <- d
  return(inst)
})


CAMERA$set("public", "plot_regional_instruments", function(region, instrument_regions=self$instrument_regions, meta_analysis_regions=self$instrument_fema_regions) {
  r <- instrument_regions[[region]]
  d <- meta_analysis_regions[[region]]
  r <- lapply(r, \(y) y %>% dplyr::mutate(z=abs(beta)/se))
  r$fema <- d
  temp <- lapply(names(r), \(y) r[[y]] %>% dplyr::mutate(pop = y)) %>% dplyr::bind_rows() %>% dplyr::select(position, z, p, pop) %>% mutate(region=region)
  th <- temp %>% dplyr::group_by(pop) %>% dplyr::arrange(desc(z)) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup()
  thf <- subset(temp, position==subset(th, pop=="fema")$position)
  ggplot2::ggplot(temp, ggplot2::aes(x=position, y=-log10(p))) +
    ggplot2::geom_point() +
    ggplot2::geom_segment(data=th, aes(x=position, xend=position, y=0, yend=-log10(p)), colour="pink") +
    ggplot2::geom_point(data=th, colour="red") +
    ggplot2::geom_segment(data=thf, aes(x=position, xend=position, y=0, yend=-log10(p)), colour="grey") +
    ggplot2::geom_point(data=thf, colour="blue") +
    ggplot2::facet_grid(pop ~ region)
})

