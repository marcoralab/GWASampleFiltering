#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

old <- commandArgs(TRUE)[1]
new <- commandArgs(TRUE)[2]
out <- commandArgs(TRUE)[3]

col.n <- c("FID", "IID", "PID", "MID", "Sex", "Phe")
col.nn <- c("none", "FIDIID", "PID_new", "MID_new", "Sex_new", "Phe_new")
col.t <- "ccccii"

new_fam <- read_table2(new, col_names = col.nn, col_types = col.t) %>%
  separate(FIDIID, c("FID_new", "IID_new"), sep = "_",
           remove = F, extra = "merge")

read_table2(old, col_names = col.n, col_types = col.t) %>%
  unite("FIDIID", FID, IID, sep = "_", remove = F) %>%
  right_join(new_fam, by = "FIDIID") %>%
  mutate(IID = ifelse(is.na(IID), IID_new, IID)) %>%
  mutate(FID = ifelse(is.na(FID), FID_new, FID)) %>%
  mutate(PID = ifelse(is.na(PID), PID_new, PID)) %>%
  mutate(MID = ifelse(is.na(MID), MID_new, MID)) %>%
  mutate(Sex = ifelse(is.na(Sex), Sex_new, Sex)) %>%
  mutate(Phe = ifelse(is.na(Phe), Phe_new, Phe)) %>%
  select(!!col.n) %>%
  write_tsv(out, col_names = F)
