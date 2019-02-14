#!/usr/bin/env Rscript

require(tibble)
require(readr)
require(dplyr)

kingstem <- commandArgs(T)[1]
exclude <- commandArgs(T)[2]

excluded <- read_tsv(exclude, col_names = c("FID", "IID"), col_types = "cc")

if ( nrow(excluded) == 0 ) {
  system(sprintf("ln -s %s.kin0 %s.popfilt.kin0", kingstem, kingstem))
  system(sprintf("ln -s %s.kin %s.popfilt.kin", kingstem, kingstem))
  q(save = "no", status = 0)
}
kincols <- cols(

  .default = col_double(),
  FID = col_character(),
  ID1 = col_character(),
  ID2 = col_character(),
  N_SNP = col_integer(),
  InfType = col_character()
)

read_tsv(paste0(kingstem, ".kin"), col_types = kincols) %>%
  anti_join(excluded, by = c("FID", "ID1" = "IID")) %>%
  anti_join(excluded, by = c("FID", "ID2" = "IID")) %>%
  write_tsv(paste0(kingstem, ".popfilt.kin"))

kin0cols <- cols(
  .default = col_double(),
  FID1 = col_character(),
  ID1 = col_character(),
  FID2 = col_character(),
  ID2 = col_character(),
  N_SNP = col_integer(),
  InfType = col_character()
)

read_tsv(paste0(kingstem, ".kin0"), col_types = kin0cols) %>%
  anti_join(excluded, by = c("FID1" = "FID", "ID1" = "IID")) %>%
  anti_join(excluded, by = c("FID2" = "FID", "ID2" = "IID")) %>%
  write_tsv(paste0(kingstem, ".popfilt.kin0"))

