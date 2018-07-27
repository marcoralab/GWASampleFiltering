#!/usr/bin/env Rscript

library(dplyr)
library(readr)

col.n = c("FID", "IID", "PID", "MID", "Sex", "Phe")
col.t = "ccccii"

fam1 <- read_table2(commandArgs(TRUE)[1],
                    col_names = col.n, col_types = col.t) %>%
  unite(FIDIID, FID, IID, sep = '_', remove = F) %>%
  rename(OLD.I = IID, OLD.F = FID)
newfam <- read_table2(commandArgs(TRUE)[2],
                      col_names = col.n, col_types = col.t) %>%
  unite(FIDIID, FID, IID, sep = '_', remove = F)

newfam %>%
  left_join(fam1, by = FIDIID) %>%
  mutate(IID = ifelse(is.na(OLD.I), IID, OLD.I)) %>%
  mutate(FID = ifelse(is.na(OLD.F), FID, OLD.F)) %>%
  write_tsv(commandArgs(TRUE)[3], col_names = F)
