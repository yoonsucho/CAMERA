
#' Fixed effects meta analysis vectorised across multiple SNPs
#' 
#' Assumes effects across studies are all on the same scale
#' 
#' @param beta_mat Matrix of betas - rows are SNPs, columns are studies
#' @param se_mat Matrix of SEs - rows are SNPs, columns are studies
#' 
#' @return list of meta analysis betas and SEs
#' @export
fixed_effects_meta_analysis <- function(beta_mat, se_mat) {
  w <- 1/se_mat^2
  beta <- rowSums(beta_mat * w) / rowSums(w)
  se <- sqrt(1/rowSums(w))
  pval <- pnorm(abs(beta/se), lower.tail=FALSE)
  return(pval)
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
  z_mat <- beta_mat / se_mat
  w_mat <- t(t(sqrt(eaf_mat * (1-eaf_mat) * 2)) * sqrt(n))
  zw <- rowSums(z_mat * w_mat) / sqrt(rowSums(w_mat^2))
  pval <- qnorm(abs(zw), lower.tail=FALSE)
  return(pval)
}



fema_regional_instruments = function(instrument_regions=self$instrument_regions, instrument_raw=self$instrument_raw)
  {
    # remove regions that have no data
    rem <- lengths(instrument_regions) == 0
    if(any(rem))
    {
      message(sum(rem), " regions have no data")
      instrument_regions <- instrument_regions[!rem]
    }
    

  }



