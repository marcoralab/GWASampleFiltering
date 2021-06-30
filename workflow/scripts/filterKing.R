#!/usr/bin/env Rscript

require(tibble)
require(readr)
require(dplyr)

kingstem <- snakemake@params[['indat']]
exclude <- snakemake@input[['exclude']]
outfile <- snakemake@output[[0]]

excluded <- read_tsv(exclude, col_names = c("FID", "IID"), col_types = "cc")

kin0file <- paste0(kingstem, ".kin0")
kinfile <- paste0(kingstem, ".kin")

if (file.exists(kinfile)) {
  message(sprintf("Loading %s", kinfile))
  if ( nrow(excluded) == 0 ) {
    system(sprintf("ln -s %s.kin %s.popfilt.kin", kingstem, kingstem))
    q(save = "no", status = 0)
  }
  kincols <- cols(
    .default = col_double(),
    FID = col_character(),
    ID1 = col_character(),
    ID2 = col_character(),
    N_SNP = col_integer()
  )

  read_tsv(kinfile, col_types = kincols) %>%
    anti_join(excluded, by = c("FID", "ID1" = "IID")) %>%
    anti_join(excluded, by = c("FID", "ID2" = "IID")) %>%
    write_tsv(paste0(kingstem, ".popfilt.kin"))
  kinout <- paste0(kingstem, ".popfilt.kin")
} else {
  message(sprintf("%s does not exist.", kin0file))
  kinout <- ""
}

if (file.exists(kin0file)) {
  message(sprintf("Loading %s", kin0file))
  if ( nrow(excluded) == 0 ) {
    system(sprintf("ln -s %s.kin0 %s.popfilt.kin0", kingstem, kingstem))
    q(save = "no", status = 0)
  }
  kin0cols <- cols(
    .default = col_double(),
    FID1 = col_character(),
    ID1 = col_character(),
    FID2 = col_character(),
    ID2 = col_character(),
    N_SNP = col_integer()
  )

  read_tsv(kin0file, col_types = kin0cols) %>%
    anti_join(excluded, by = c("FID1" = "FID", "ID1" = "IID")) %>%
    anti_join(excluded, by = c("FID2" = "FID", "ID2" = "IID")) %>%
    write_tsv(paste0(kingstem, ".popfilt.kin0"))
  kin0out <- paste0(kingstem, ".popfilt.kin0")
} else {
  message(sprintf("%s does not exist.", kin0file))
  kin0out <- ""
}

fileConn <- file(outfile)
writeLines(c(kinout,kin0out), fileConn)
close(fileConn)
