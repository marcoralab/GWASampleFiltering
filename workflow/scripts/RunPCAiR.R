#!/usr/bin/env Rscript

require(GENESIS)
library(SNPRelate)
require(tibble)
require(readr)
require(dplyr)

fstem <- "output/ADGC/x_present_AA/REAAADI-AA_filtered_PCApre"
kingstem <-"output/ADGC/x_present_AA/REAAADI-AA_IBDQC.all.popfilt"
iterlist <- "output/ADGC/x_present_AA/REAAADI-AA_exclude.relatedness"

fstem <- snakemake@params[['stem']]
kingstem <- snakemake@params[['king']]
iterlist <- snakemake@input[['iterative']]

fam <- read_table2(paste0(fstem, ".fam"), col_types = "ccccii",
                   col_names = c("FID", "IID", "PID", "MID", "Sex", "Pheno"))

iterlist <- iterlist %>%
  read_tsv(col_types = "cc")

iids <- fam$IID

use_iterlist <- function(fam, iterlist) {
  if ( nrow(iterlist) == 0 ) {
    unrel <- fam %>%
      select(FID, IID)
  } else {
    unrel <- fam %>%
      anti_join(iterlist, by = c("FID", "IID")) %>%
      select(FID, IID)
  }
  rel <- iterlist
  return(list(rels = rel, unrels = unrel))
}

import_KING_external <- function(fam, kingstem) {
  if (length(unique(fam$FID)) == 1) {
    kingext <- ".kin"
  } else if (length(unique(fam$FID)) == length(fam$FID)) {
    kingext <- ".kin0"
  } else {
    kingext <- c(".kin0", ".kin")
  }
  kingfiles <- paste0(kingstem, kingext)
  if (exists("kingToMatrix")) {
    KINGmat <- kingToMatrix(kingfiles, sample.include = fam$IID, estimator = "Kinship")
  } else if (length(kingfiles) == 2) {
    KINGmat <- king2mat(file.kin0 = kingfiles[1],
                        file.kin = kingfiles[2],
                        iids = fam$IID)
  } else {
    if (length(kingfiles) == 1 & !grepl("kin0", kingfiles)) {
      warning("This version of GENESIS requires a kin0 file, but everyone has the same FID. Generating a kin0 file with the same pairs as in the kin file.")
      read_table2(kingfiles, col_types = cols(
        .default = col_double(),
        FID = col_character(),
        ID1 = col_character(),
        ID2 = col_character(),
        N_SNP = col_integer(),
        InfType = col_character()
      )) %>%
        mutate(FID2 = FID) %>%
        rename(FID1 = FID) %>%
        select(FID1, ID1, FID2, ID2, everything()) %>%
        write_tsv(paste0(kingfiles,"0"), col_names = T)
      kingfiles <- paste0(kingfiles,"0")
    }
    KINGmat <- king2mat(file.kin0 = kingfiles, iids = fam$IID)
  }
  return(KINGmat)
}

pcair_segment <- function(KINGmat, divthresh="default") {
  if (divthresh == "default") {
    if (exists("kingToMatrix")) {
      raws <- pcairPartition(kinobj = KINGmat, divobj = KINGmat)
    } else {
      raws <- pcairPartition(kinMat = KINGmat, divMat = KINGmat)
    }
  } else {
    warning("Initial div threshold failed to find unrelated set.")
    message(sprintf("Retrying with threshold of %s.", divthresh))
    if (exists("kingToMatrix")) {
      raws <- pcairPartition(kinobj = KINGmat, divobj = KINGmat, div.thresh = divthresh)
    } else {
      raws <- pcairPartition(kinMat = KINGmat, divMat = KINGmat, div.thresh = divthresh)
    }
  }
  return(raws)
}

tryDivthresh <- function(KINGmat, thresholds, retry=F) {
  ret <- tryCatch({
    if (retry) {
      raws <- pcair_segment(KINGmat, thresholds[1])
    } else {
      raws <- pcair_segment(KINGmat)
    }
    raws
  }, error = function (e) {
    if ( grepl("must be the same length as the vector", e) && length(thresholds) > 1 ) {
      ret <- tryDivthresh(KINGmat, thresholds[2:length(thresholds)], retry=T)
      return(ret)
    } else {
      stop(e)
    }
  })
  return(ret)
}

use_pcair <- function(fam, kingstem) {
  raws <- import_KING_external(fam, kingstem) %>%
    tryDivthresh(-2^-c(6, 6.5, 7, 8))
  print("Related set:")
  print(raws$rels)
  print("Unrelated set:")
  print(raws$unrels)
  if ( length(raws$rels) + length(raws$unrels) < nrow(fam) ) {
    stop("Length mismatch between KING output and famfile.")
  }
  FIDIID <- fam %>%
    select(FID, IID)
  if (is.null(raws$rels)) raws$rels <- c("")
  rel <- raws$rels %>%
    tibble(IID = .) %>%
    inner_join(FIDIID, by = "IID") %>%
    select(FID, IID)
  unrel <- raws$unrels %>%
    tibble(IID = .) %>%
    inner_join(FIDIID, by = "IID") %>%
    select(FID, IID)
  return(list(rels = rel, unrels = unrel))
}

pca_intern <- function(gds, unrel_iid) {
  PCA_unrelate <- snpgdsPCA(gds, sample.id = unrel_iid)
  
  loadings <- snpgdsPCASNPLoading(PCA_unrelate, gds)
  
  PCA_all <- snpgdsPCASampLoading(loadings, gds)
  
  PCA_tibble <- function(pc, iid) {
    pc = as.data.frame(pc)
    colnames(pc) <- paste0("PC", 1:ncol(pc))
    as_tibble(pc) %>%
      mutate(IID = as.character(iid)) %>%
      select(IID, everything())
  }
  
  PCA_all_fmt <- PCA_tibble(PCA_all$eigenvect, PCA_all$sample.id)
  PCA_eigenvals <- PCA_unrelate$eigenval[1:ncol(PCA_unrelate$eigenvect)]
  
  return(list(eigenvals = PCA_eigenvals,
              eigenvecs = PCA_all_fmt,
              PCA_unrel = PCA_unrelate,
              PCA_all   = PCA_all))
}

partition <- function(fam, kingstem, iterlist, iids) {
  if ( length(unique(iids)) < length(iids)) {
    parts <- use_iterlist(fam, iterlist)
    method <- "iterative (plink)"
    reason <- "There are non-unique IIDs."
    out <- list(parts = parts, method = method, reason = reason,
                iid_repeats = TRUE)
  } else {
    out <- tryCatch({
      parts <- use_pcair(fam, kingstem)
      method <- "PC-AiR"
      reason <- "This is the default method."
      list(parts = parts, method = method, reason = reason, iid_repeats = FALSE)
    }, error = function(e) {
      if (grepl("unmatched node provided", e)) {
        parts <- use_iterlist(fam, iterlist)
        method <- "iterative (SNPRelate)"
        reason <- "PC-AiR was unable to partition based on relatedness."
        list(parts = parts, method = method, reason = reason,
             iid_repeats = FALSE)
      } else {
        stop(e)
      }
    })
  }
  return(out)
}

partitioning <- partition(fam, kingstem, iterlist, iids)
parts <- partitioning$parts
method <- partitioning$method
reason <- partitioning$reason
iid_repeats <- partitioning$iid_repeats
unrel <- parts$unrels %>%
  mutate(cluster = "unrelated")
rel <- parts$rels %>%
  mutate(cluster = "related")
relatedness <- bind_rows(unrel, rel)
logtxt <- c(method, reason)

write_delim(relatedness, paste0(fstem, ".unrel"), col_names = F)
write_lines(logtxt, paste0(fstem, ".partition.log"))

if (iid_repeats) {
  save(iid_repeats, method, reason, relatedness, file = paste0(fstem, ".rda"))
} else {
  snpgdsBED2GDS(
    paste0(fstem, ".bed"), paste0(fstem, ".fam"), paste0(fstem, ".bim"),
    paste0(fstem, ".gds"), family=T)
  gds <- snpgdsOpen(paste0(fstem, ".gds"))
  pca <- pca_intern(gds, parts$unrels$IID)
  pca$eigenvecs <- pca$eigenvecs %>%
    left_join(select(fam, FID, IID), by = "IID") %>%
    select(FID, IID, everything())
  snpgdsClose(gds)
  eigenvals <- pca$eigenvals
  eigenvecs <- pca$eigenvecs
  write_lines(eigenvals, paste0(fstem, ".pcair.eigenval"))
  write_delim(eigenvecs, paste0(fstem, ".pcair.eigenvec"), col_names = F)
  save(iid_repeats, method, reason, relatedness, eigenvals, eigenvecs,
       file = paste0(fstem, ".rda"))
}
