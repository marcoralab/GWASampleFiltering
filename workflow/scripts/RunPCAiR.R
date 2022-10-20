#!/usr/bin/env Rscript

require(GENESIS)
library(SNPRelate)
require(tibble)
require(readr)
suppressMessages(require(dplyr))

outdir <- "results/intermediate/post-impute_filter/output"
fstem <- paste0(outdir, "/imputed_filtered_PCApre")
fstem_unpruned <- paste0(outdir, "/imputed_callRate")
kingstem <- paste0(outdir, "/imputed_IBDQC.all.popfilt")
iterlist <- paste0(outdir, "/imputed_exclude.relatedness")

fstem <- snakemake@params[["stem"]]
fstem_unpruned <- snakemake@params[["stem_unpruned"]]
kingstem <- snakemake@params[["king"]]
iterlist <- snakemake@input[["iterative"]]

fam <- read_table2(paste0(fstem, ".fam"), col_types = "ccccii",
                   col_names = c("FID", "IID", "PID", "MID", "Sex", "Pheno"))

iterlist <- iterlist %>%
  read_tsv(col_types = "cc")

iids <- fam$IID

use_iterlist <- function(fam, iterlist) {
  if (nrow(iterlist) == 0) {
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

import_king_external <- function(fam, kingstem) {
  if (length(unique(fam$FID)) == 1) {
    kingext <- ".kin"
  } else if (length(unique(fam$FID)) == length(fam$FID)) {
    kingext <- ".kin0"
  } else {
    kingext <- c(".kin0", ".kin")
  }
  kingfiles <- paste0(kingstem, kingext)
  if (exists("kingToMatrix")) {
    kingmat <- GENESIS::kingToMatrix(kingfiles,
                                     sample.include = fam$IID,
                                     estimator = "Kinship")
  } else if (length(kingfiles) == 2) {
    kingmat <- GENESIS::king2mat(file.kin0 = kingfiles[1],
                                 file.kin = kingfiles[2],
                                 iids = fam$IID)
  } else {
    if (length(kingfiles) == 1 && !grepl("kin0", kingfiles)) {
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
        write_tsv(paste0(kingfiles, "0"), col_names = TRUE)
      kingfiles <- paste0(kingfiles, "0")
    }
    kingmat <- GENESIS::king2mat(file.kin0 = kingfiles, iids = fam$IID)
  }
  return(kingmat)
}

generate_king_intern <- function(fam, fstem_unpruned) {
  gdsfile <- paste0(fstem_unpruned, ".gds")
  snpgdsBED2GDS(
    paste0(fstem_unpruned, ".bed"),
    paste0(fstem_unpruned, ".fam"),
    paste0(fstem_unpruned, ".bim"),
    gdsfile, family = TRUE)
  gds <- snpgdsOpen(gdsfile)
  ibd_robust <- snpgdsIBDKING(gds, num.thread = 48)
  snpgdsClose(gds)
  GENESIS::kingToMatrix(ibd_robust, sample.include = fam$IID) %>%
    return()
}

load_king_intern <- function(fam, rds) {
  read_rds(rds) %>%
    GENESIS::kingToMatrix(sample.include = fam$IID) %>%
    return()
}

pcair_segment <- function(kinmat, divmat, divthresh = "default") {
  if (divthresh == "default") {
    if (exists("kingToMatrix")) {
      raws <- GENESIS::pcairPartition(kinobj = kinmat, divobj = divmat)
    } else {
      raws <- GENESIS::pcairPartition(kinMat = kinmat, divMat = divmat)
    }
  } else {
    warning("Initial div threshold failed to find unrelated set.")
    message(sprintf("Retrying with threshold of %s.", divthresh))
    if (exists("kingToMatrix")) {
      raws <- GENESIS::pcairPartition(kinobj = kinmat, divobj = divmat,
                                      div.thresh = divthresh)
    } else {
      raws <- GENESIS::pcairPartition(kinMat = kinmat, divMat = divmat,
                                      div.thresh = divthresh)
    }
  }
  return(raws)
}

try_divthresh <- function(kinmat, divmat, thresholds, retry=F) {
  ret <- tryCatch({
    if (retry) {
      raws <- pcair_segment(kinmat, divmat, thresholds[1])
    } else {
      raws <- pcair_segment(kinmat, divmat)
    }
    raws
  }, error = function(e) {
    if (grepl("must be the same length as the vector", e) &&
        length(thresholds) > 1) {
      ret <- try_divthresh(kinmat, divmat, thresholds[2:length(thresholds)],
                          retry = TRUE)
      return(ret)
    } else {
      stop(e)
    }
  })
  return(ret)
}

use_pcair <- function(fam, kinmat, divmat) {
  raws <- try_divthresh(kinmat, divmat, -2^-c(6, 6.5, 7, 8))
  print("Related set:")
  print(raws$rels)
  print("Unrelated set:")
  print(raws$unrels)
  if (length(raws$rels) + length(raws$unrels) < nrow(fam)) {
    stop("Length mismatch between KING output and famfile.")
  }
  fidiid <- fam %>%
    select(FID, IID)
  if (is.null(raws$rels)) raws$rels <- c("")
  rel <- raws$rels %>%
    tibble(IID = .) %>%
    inner_join(fidiid, by = "IID") %>%
    select(FID, IID)
  unrel <- raws$unrels %>%
    tibble(IID = .) %>%
    inner_join(fidiid, by = "IID") %>%
    select(FID, IID)
  return(list(rels = rel, unrels = unrel))
}

pca_intern <- function(gds, unrel_iid, fast = FALSE) {
  if (fast) {
    algo <- "randomized"
  } else {
    algo <- "exact"
  }
  pca_unrelate <- snpgdsPCA(gds, sample.id = unrel_iid, num.thread = 48,
                            algorithm = algo)
  loadings <- snpgdsPCASNPLoading(pca_unrelate, gds, num.thread = 48)
  pca_all <- snpgdsPCASampLoading(loadings, gds, num.thread = 48)

  pca_tibble <- function(pc, iid) {
    pc <- as.data.frame(pc)
    colnames(pc) <- paste0("PC", seq_len(ncol(pc)))
    as_tibble(pc) %>%
      mutate(IID = as.character(iid)) %>%
      select(IID, everything())
  }

  pca_all_fmt <- pca_tibble(pca_all$eigenvect, pca_all$sample.id)
  pca_eigenvals <- pca_unrelate$eigenval[seq_len(ncol(pca_unrelate$eigenvect))]

  return(list(eigenvals = pca_eigenvals,
              eigenvecs = pca_all_fmt,
              PCA_unrel = pca_unrelate,
              PCA_all   = pca_all))
}

partition <- function(fam, kinmat, divmat, iterlist, iids) {
  if (length(unique(iids)) < length(iids)) {
    parts <- use_iterlist(fam, iterlist)
    method <- "iterative (plink)"
    reason <- "There are non-unique IIDs."
    out <- list(parts = parts, method = method, reason = reason,
                iid_repeats = TRUE)
  } else {
    out <- tryCatch({
      parts <- use_pcair(fam, kinmat, divmat)
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
  unrel <- out$parts$unrels %>%
    mutate(cluster = "unrelated")
  rel <- out$parts$rels %>%
    mutate(cluster = "related")
  out$relatedness <- bind_rows(unrel, rel)
  return(out)
}

choose_eigencount <- function(eigenvals, thresh) {
  eigenvals_no_na <- tidyr::replace_na(eigenvals, 0)
  tot <- sum(eigenvals_no_na)
  pve <- 0
  n_eig <- 0
  while (pve < thresh && n_eig <= length(eigenvals_no_na)) {
    n_eig <- n_eig + 1
    pve <- sum(eigenvals_no_na[seq_len(n_eig)] / tot)
  }
  message(sprintf("Selecting %i PCs out of %i. (PVE = %f)",
                  n_eig, length(eigenvals_no_na), pve))
  return(n_eig)
}

if (length(unique(iids)) < length(iids)) {
  kingobj <- ""
} else if (length(iids) < 20000) {
  kingobj <- import_king_external(fam, kingstem)
} else {
  kingobj <- generate_king_intern(fam, fstem_unpruned)
}

message("Performing initial partitioning")
partitioning <- partition(fam, kingobj, kingobj, iterlist, iids)

parts <- partitioning$parts
method <- partitioning$method
reason <- partitioning$reason
iid_repeats <- partitioning$iid_repeats
relatedness <- partitioning$relatedness

write_lines(c(method, reason), paste0(fstem, ".partition.log"))

save.image(file = paste0(fstem, ".rda"))

if (iid_repeats) {
  write_delim(relatedness, paste0(fstem, ".unrel"), col_names = FALSE)
  save(iid_repeats, method, reason, relatedness, file = paste0(fstem, ".rda"))
} else {
  message("Creating Pruned GDS")
  gdsfile <- paste0(fstem, ".gds")
  snpgdsBED2GDS(
    paste0(fstem, ".bed"), paste0(fstem, ".fam"), paste0(fstem, ".bim"),
    gdsfile, family = TRUE)

  parobj <- BiocParallel::MulticoreParam(workers = 48)

  # first iteration
  save.image(file = paste0(fstem, ".rda"))
  message("Performing initial PCA")
  gds <- snpgdsOpen(gdsfile)
  pca_initial <- pca_intern(gds, parts$unrels$IID, fast = TRUE)
  n_eig <- choose_eigencount(pca_initial$PCA_unrel$eigenval, 0.8)
  pcmat_initial <- pca_initial$PCA_all$eigenvect
  rownames(pcmat_initial) <- as.character(pca_initial$PCA_all$sample.id)
  snpgdsClose(gds)
  # create a GenotypeData class object
  save.image(file = paste0(fstem, ".rda"))

  message("Performing initial PCRelate")
  gds_geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)
  genodata <- GWASTools::GenotypeData(gds_geno)
  genoiter <- GWASTools::GenotypeBlockIterator(genodata)
  pcrel_initial <- pcrelate(gdsobj = genoiter,
                            pcs = pcmat_initial[, seq_len(min(n_eig, 16))],
                            training.set = parts$unrels$IID,
                            small.samp.correct = nrow(fam) < 5000,
                            ibd.probs = FALSE,
                            BPPARAM = parobj)
  GWASTools::close(genodata)
  # second iteration

  message("Repartitioning")
  kinmat <- pcrelateToMatrix(pcrel_initial, thresh = 2^(-11 / 2))
  partitioning_second <- partition(fam, kinmat, kingobj, iterlist, iids)

  message("Performing final PCA")
  gds <- snpgdsOpen(gdsfile)
  pca <- pca_intern(gds, partitioning_second$parts$unrels$IID)
  snpgdsClose(gds)

  pcmat <- pca$PCA_all$eigenvect
  rownames(pcmat) <- as.character(pca$PCA_all$sample.id)
  n_eig_final <- choose_eigencount(pca$PCA_unrel$eigenval, 0.8)

  message("Performing final PCRelate")
  gds_geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)
  genodata <- GWASTools::GenotypeData(gds_geno)
  genoiter <- GWASTools::GenotypeBlockIterator(genodata)
  pcrel <- pcrelate(gdsobj = genoiter,
                    pcs = pcmat[, seq_len(min(n_eig, 16))],
                    sample.include = iids,
                    training.set = partitioning_second$parts$unrels$IID,
                    small.samp.correct = nrow(fam) < 5000,
                    BPPARAM = parobj)
  GWASTools::close(genodata)

  message("Writing Results")

  pca$eigenvecs <- pca$eigenvecs %>%
    left_join(select(fam, FID, IID), by = "IID") %>%
    select(FID, IID, everything())
  eigenvals <- pca$eigenvals
  eigenvecs <- pca$eigenvecs
  write_lines(eigenvals, paste0(fstem, ".pcair.eigenval"))
  write_delim(eigenvecs, paste0(fstem, ".pcair.eigenvec"), col_names = FALSE)

  relatedness <- partitioning_second$relatedness
  write_delim(relatedness, paste0(fstem, ".unrel"), col_names = FALSE)

  save(iid_repeats, partitioning_second, method, reason, relatedness,
       eigenvals, eigenvecs, pcrel, file = paste0(fstem, ".rda"))
}
