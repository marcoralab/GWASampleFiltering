#!/usr/bin/env Rscript

require(GENESIS)
require(tibble)
require(readr)
require(dplyr)

fstem <- commandArgs(T)[1]
kingstem <- commandArgs(T)[2]
iterlist <- commandArgs(T)[3]

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

use_pcair <- function(fam, kingstem) {
  KINGmat <- king2mat(file.kin0 = paste0(kingstem, ".kin0"),
                      file.kin = paste0(kingstem, ".kin"),
                      iids = fam$IID)
  raws <- pcairPartition(kinMat = KINGmat, divMat = KINGmat)
  if ( length(raws$rels) + length(raws$unrels) < nrow(fam) ) {
    stop("Length mismatch between KING output and famfile.")
  }
  FIDIID <- fam %>%
    select(FID, IID)
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

