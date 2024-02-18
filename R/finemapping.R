# for each instrument region + ld matrix we can perform susie finemapping
# do this independently in each population
# find credible sets that overlap - to try to determine the best SNP in a region to be used as instrument
#' @description
#' Fine-mapping using susieR (https://github.com/stephenslab/susieR) to extract instruments for the exposure for multiple populations. The function identifies credible sets that overlap between the populations and determine the best SNP in a genomic region to be used as instrument.
#' @param dat Genomic regions identified by using \code{x$extract_instrument_regions()}
#' @param ld LD matrix obtained by using \code{x$regional_ld_matrices()}
#' @return Result from susieR in x$susie_results. Data frame in x$instrument_susie
CAMERA$set("public", "susie_finemap_regions", function(dat = self$instrument_regions, ld = self$ld_matrices) {
  resnps <- names(ld)
  exp <- self$exposure_ids
  susie <- lapply(1:length(resnps), function(r) {
    tryCatch(
      {
        snp <- rownames(ld[[r]][[exp[1]]])[rownames(ld[[r]][[exp[1]]]) %in% rownames(ld[[r]][[exp[2]]])]
        snp <- strsplit(snp, "_") %>% sapply(., function(x) x[1])

        d <- lapply(1:length(exp), function(i) {
          index <- which(dat[[r]][[i]]$rsid %in% snp)
          d <- dat[[r]][[i]][index, ]
          return(d)
        })

        su <- lapply(1:length(exp), function(i) {
          message("Running susie for ", exp[i])
          susieR::susie_rss(z = d[[i]]$beta / d[[i]]$se, R = ld[[r]][[i]])
        })

        message("Finding overlaps [", r, "/", length(resnps), "]")
        suo <- private$susie_overlaps(su[[1]], su[[2]])

        out <- list(
          chr = d[[1]]$chr,
          position = d[[1]]$position,
          radius = self$radius,
          a1 = d[[1]],
          a2 = d[[2]],
          su1 = su[[1]],
          su2 = su[[2]],
          suo = suo
        )

        null <- c(is.null(su[[1]]$sets$cs), is.null(su[[2]]$sets$cs))
        if (null[1] & !null[2]) {
          out$type <- "pop2"
        } else if (null[2] & !null[1]) {
          out$type <- "pop1"
        } else if (!null[1] & !null[2]) {
          out$type <- "shared"
        } else {
          out$type <- "drop"
        }

        if (out$type == "shared") {
          if (nrow(suo) == 0) {
            out$cs_overlap <- FALSE
            temp <- dplyr::inner_join(d[[1]], d[[2]], by = "rsid") %>%
              dplyr::mutate(pvalrank = rank(p.x) / nrow(d[[1]]) + rank(p.y) / nrow(d[[1]])) %>%
              dplyr::arrange(pvalrank)
            out$bestsnp <- temp$rsid[1]
          } else {
            out$cs_overlap <- TRUE
            out$bestsnp <- d[[1]]$rsid[suo$v[1]]
          }
        }

        if (out$type == "pop1") {
          out$cs_overlap <- FALSE
          out$bestsnp <- d[[1]] %>%
            dplyr::arrange(p) %>%
            {
              .$rsid[1]
            }
        }
        if (out$type == "pop2") {
          out$cs_overlap <- FALSE
          out$bestsnp <- d[[2]] %>%
            dplyr::arrange(p) %>%
            {
              .$rsid[1]
            }
        }
        return(out)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })

  names(susie) <- names(self$instrument_regions)
  self$susie_results <- susie

  # susie <- susie[!sapply(susie, is.null)]
  o <- unique(lapply(resnps, function(r) {
    susie[[r]]$bestsnp
  }) %>% unlist())
  instrument_susie <- lapply(resnps, function(r) {
    tryCatch(
      {
        lapply(exp, function(id) {
          dat[[r]][[id]] %>%
            subset(., rsid %in% o) %>%
            dplyr::bind_rows() %>%
            dplyr::arrange(id, chr, position)
        })
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()

  t <- self$instrument_raw %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::select(id, units, samplesize)

  instrument_susie <- dplyr::left_join(instrument_susie, t, by = "id") %>% as.data.frame()
  instrument_susie <- lapply(exp, function(i) {
    subset(instrument_susie, id == i) %>% dplyr::distinct(., rsid, .keep_all = TRUE)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(method = "susie")
  self$instrument_susie <- instrument_susie
  invisible(self)
})


#' @importFrom dplyr bind_rows mutate arrange
#' @importFrom tibble tibble
CAMERA$set("private", "susie_overlaps", function(su1, su2) {
  l <- list()
  k <- 1
  s1 <- su1$sets$cs
  s2 <- su2$sets$cs
  for (i in 1:length(s1))
  {
    for (j in 1:length(s2))
    {
      if (any(s1[[i]] %in% s2[[j]])) {
        ind <- s1[[i]] %in% s2[[j]]
        v <- s1[[i]][ind]
        l[[k]] <- tibble::tibble(s1 = i, s2 = j, v = v)
        k <- k + 1
      }
    }
  }
  l <- dplyr::bind_rows(l)
  if (nrow(l) > 0) {
    l$pip1 <- su1$pip[l$v]
    l$pip2 <- su2$pip[l$v]
    l$piprank1 <- rank(-su1$pip)[l$v] / length(su1$pip)
    l$piprank2 <- rank(-su2$pip)[l$v] / length(su2$pip)
    l <- l %>%
      dplyr::mutate(piprank = piprank1 + piprank2) %>%
      dplyr::arrange(piprank)
  }
  return(l)
})

# PAINTOR allows finemapping jointly across multiple populations
# returns a posterior probability of inclusion for each SNP
#' @description
#' Fine-mapping using PAINTOR (https://github.com/gkichaev/PAINTOR_V3.0) to extract instruments for the exposure for multiple populations. The function chooses the SNP with highest posterior probability and associations in each exposure.
#' @param region Genomic regions identified by using \code{x$extract_instrument_regions()}
#' @param ld LD matrix obtained by using \code{x$regional_ld_matrices()}
#' @param PAINTOR Path to executable PAINTOR. Default="PAINTOR"
#' @param workdir Working directory to save the output. Default=tempdir()
#' @return Result from PAINTOR in x$paintor_results. Data frame in x$instrument_paintor
CAMERA$set("public", "paintor_finemap_regions", function(region = self$instrument_regions, ld = self$ld_matrices, PAINTOR = "PAINTOR", workdir = tempdir()) {
  id <- self$exposure_ids
  nid <- length(region)
  snps <- lapply(1:nid, function(i) {
    region[[i]][[1]]$rsid
  })
  zs <- lapply(1:nid, function(i) {
    o <- list()
    lapply(1:length(id), function(id) {
      o[[paste0("ZSCORE.P", id)]] <- region[[i]][[id]]$beta / region[[i]][[id]]$se
      return(tibble::as_tibble(o))
    }) %>% dplyr::bind_cols()
  })

  locus <- lapply(1:nid, function(i) {
    tryCatch(
      {
        l <- list()
        l <- tibble::tibble(CHR = region[[i]][[1]]$chr, POS = region[[i]][[1]]$position, RSID = region[[i]][[1]]$rsid)
        l <- l %>% dplyr::bind_cols(., zs[[i]])
        ldsnp <- strsplit(rownames(ld[[i]][[1]]), "_") %>% sapply(., function(x) x[1])
        snp <- l$RSID %in% ldsnp
        index <- which(l$RSID %in% ldsnp)
        d <- l[index, ]
        return(d)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })

  anno <- lapply(1:nid, function(i) {
    tryCatch(
      {
        tibble::tibble(null = rep(1, nrow(locus[[i]])))
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })

  # write.files
  lapply(1:nid, function(i) {
    tryCatch(
      {
        write.table(locus[i][[1]], file = file.path(workdir, paste0("Locus", i)), row = F, col = T, qu = F)
        write.table(anno[[i]], file = file.path(workdir, paste0("Locus", i, ".annotations")), row = F, col = T, qu = F)
        write.table(paste0("Locus", i), file = file.path(workdir, paste0("input.files", i)), row = F, col = F, qu = F)

        lapply(1:length(id), function(id) {
          write.table(ld[[i]][[id]], file = file.path(workdir, paste0("Locus", i, ".LD", id)), row = FALSE, col = FALSE, qu = FALSE)
        })
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })

  LDname <- paste0("LD", 1:length(id), collapse = ",")
  Zhead <- paste(names(zs[[1]]), collapse = ",")

  wd <- getwd()
  setwd(workdir)
  lapply(1:nid, function(i) {
    tryCatch(
      {
        system(glue::glue("{PAINTOR} -input input.files{i} -Zhead {Zhead} -LDname {LDname} -in {workdir}/ -out {workdir}/ -mcmc -annotations null"))
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })

  res <- lapply(1:nid, function(i) {
    tryCatch(
      {
        data.table::fread(file.path(workdir, paste0("Locus", i, ".results"))) %>%
          dplyr::mutate(ZSCORE.SUM = ZSCORE.P1 + ZSCORE.P2) %>%
          dplyr::rename(., chr = CHR, position = POS, rsid = RSID)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })
  names(res) <- names(self$instrument_regions)
  unlink(workdir)
  setwd(wd)
  self$paintor_results <- res

  # res <- res[!sapply(res, is.null)]
  bestsnp <- lapply(1:length(res), function(r) {
    tryCatch(
      {
        if (sum(res[[r]]$Posterior_Prob) == 0) {
          NULL
        } else {
          res[[r]] %>%
            dplyr::arrange(desc(Posterior_Prob)) %>%
            {
              .$rsid[1]
            }
        }
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })

  instrument_paintor <- lapply(1:nid, function(r) {
    tryCatch(
      {
        lapply(1:length(id), function(id) {
          region[[r]][[id]] %>%
            subset(., rsid %in% bestsnp) %>%
            dplyr::bind_rows() %>%
            dplyr::arrange(id, chr, position)
        })
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()

  t <- self$instrument_raw %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::select(id, units, samplesize)
  instrument_paintor <- dplyr::left_join(instrument_paintor, t, by = "id") %>% as.data.frame()
  instrument_paintor <- lapply(id, function(i) {
    tryCatch(
      {
        subset(instrument_paintor, id == i) %>% dplyr::distinct(., rsid, .keep_all = TRUE)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(method = "paintor")
  self$instrument_paintor <- instrument_paintor
  invisible(self)
})


#' @description
#' Fine-mapping using MsCAVIAR (https://github.com/nlapier2/MsCAVIAR) to extract instruments for the exposure for multiple populations.
#' @param region Genomic regions identified by using \code{x$extract_instrument_regions()}
#' @param ld LD matrix obtained by using \code{x$regional_ld_matrices()}
#' @param MsCAVIAR Path to executable MsCAVIAR. Default="MsCAVIAR"
#' @param workdir Working directory to save the output. Default=tempdir()
#' @return Result from MsCAVIAR in x$mscaviar_results. Data frame in x$instrument_mscaviar
CAMERA$set("public", "MsCAVIAR_finemap_regions", function(region = self$instrument_regions, ld = self$ld_matrices, MsCAVIAR = "MsCAVIAR", workdir = tempdir()) {
  id <- self$exposure_ids
  nid <- length(region)
  snps <- lapply(1:nid, function(i) {
    region[[i]][[1]]$rsid
  })
  n <- paste0(unique(self$instrument_raw$samplesize), collapse = ",")

  zs <- lapply(1:nid, function(i) {
    lapply(1:length(id), function(id) {
      zs <- region[[i]][[id]] %>%
        dplyr::mutate(z = beta / se) %>%
        dplyr::select(rsid, z)
      names(zs)[2] <- c(paste0("ZSCORE.P", id))
      ldsnp <- strsplit(rownames(ld[[i]][[1]]), "_") %>% sapply(., function(x) x[1])
      snp <- zs$rsid %in% ldsnp
      index <- which(zs$rsid %in% ldsnp)
      zs <- zs[index, ]
      write.table(zs, file = file.path(workdir, paste0("z_", i, "_", id, ".zscores")), row = F, col = F, qu = F)
      return(tibble::as_tibble(zs))
    })
  })

  lapply(1:nid, function(i) {
    lapply(1:length(id), function(id) {
      write.table(ld[[i]][[id]], file = file.path(workdir, paste0("ld_", i, "_", id, ".ld")), row = FALSE, col = FALSE, qu = FALSE)
    })
  })

  thisdir <- getwd()
  setwd(workdir)

  lapply(1:nid, function(i) {
    writeLines(c(paste0("z_", i, "_", 1:length(id), ".zscores")), file(file.path(workdir, paste0("zfiles", i, ".txt"))))
    writeLines(c(paste0("ld_", i, "_", 1:length(id), ".ld")), file(file.path(workdir, paste0("ldfiles", i, ".txt"))))
  })

  lapply(1:nid, function(i) {
    system(glue::glue("{MsCAVIAR} -l ldfiles{i}.txt -z zfiles{i}.txt -n {n} -o results{i}"))
    message("Analysis [", i, "/", nid, "]")
  })

  res <- lapply(1:nid, function(i) {
    tryCatch(
      {
        res <- data.table::fread(file.path(workdir, paste0("results", i, "_post.txt")))
        names(res)[1] <- c("rsid")
        names(res)[2] <- c("Posterior_Prob")
        return(res)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  })


  zs2 <- suppressMessages(lapply(1:nid, function(i) {
    zs[[i]] %>%
      dplyr::bind_cols() %>%
      dplyr::select(rsid = rsid...1, ZSCORE.P1, ZSCORE.P2)
  }))

  res <- suppressMessages(lapply(1:nid, function(i) {
    tryCatch(
      {
        res <- res[[i]] %>%
          dplyr::left_join(., zs2[[i]]) %>%
          dplyr::left_join(., region[[i]][[1]]) %>%
          dplyr::select(rsid, chr, position, ZSCORE.P1, ZSCORE.P2, Posterior_Prob)
        return(res)
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }))
  names(res) <- names(self$instrument_regions)
  unlink(workdir)
  setwd(thisdir)

  self$mscaviar_results <- res

  bestsnp <- lapply(1:length(res), function(r) {
    if (sum(res[[r]]$Posterior_Prob) == 0) {
      NULL
    } else {
      res[[r]] %>%
        dplyr::arrange(desc(Posterior_Prob)) %>%
        {
          .$rsid[1]
        }
    }
  })

  instrument_mscaviar <- lapply(1:nid, function(r) {
    lapply(1:length(id), function(id) {
      region[[r]][[id]] %>%
        subset(., rsid %in% bestsnp) %>%
        dplyr::bind_rows() %>%
        dplyr::arrange(id, chr, position)
    })
  }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()

  t <- self$instrument_raw %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::select(id, units, samplesize)
  instrument_mscaviar <- dplyr::left_join(instrument_mscaviar, t, by = "id") %>% as.data.frame()
  instrument_mscaviar <- lapply(id, function(i) {
    subset(instrument_mscaviar, id == i)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::relocate(rsid, chr, position, id, beta, se, p, ea, nea, units, samplesize) %>%
    dplyr::mutate(method = "mscaviar")
  self$instrument_mscaviar <- instrument_mscaviar
  invisible(self)
})
