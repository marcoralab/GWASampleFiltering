#!/usr/bin/env Rscript

require(GENESIS)
require(tibble)
require(readr)
require(dplyr)

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

if ( length(unique(iids)) < length(iids) ) {
  parts <- use_iterlist(fam, iterlist)
  method <- "iterative"
  reason <- "There are non-unique IIDs."
} else {
  out <- tryCatch({
    parts <- use_pcair(fam, kingstem)
    method <- "PC-AiR"
    reason <- "This is the default method."
    list(parts = parts, method = method, reason = reason)
  }, error = function (e) {
    if ( grepl("unmatched node provided", e) ) {
      parts <- use_iterlist(fam, iterlist)
      method <- "iterative"
      reason <- "PC-AiR was unable to partition based on relatedness."
      list(parts = parts, method = method, reason = reason)
    } else {
      stop(e)
    }
  })
  parts <- out$parts
  method <- out$method
  reason <- out$reason
}

output <- parts$unrels %>%
  mutate(cluster = "unrelated")
logtxt <- c(method, reason)

write_delim(output, paste0(fstem, ".unrel"), col_names = F)
write_lines(logtxt, paste0(fstem, ".partition.log"))
