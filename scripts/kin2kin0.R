#!/usr/bin/env Rscript

suppressMessages(require(readr))
suppressMessages(require(dplyr))

genome.file <- commandArgs(TRUE)[1]

read_table2(genome.file, col_types = cols(
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
  write_tsv(paste0(genome.file,"0"), col_names = T)

