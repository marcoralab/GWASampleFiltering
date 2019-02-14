#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

old <- commandArgs(TRUE)[1]
new <- commandArgs(TRUE)[2]
out <- commandArgs(TRUE)[3]
ped <- commandArgs(TRUE)[4]

col.n <- c("FID", "IID", "PID", "MID", "Sex", "Phe")
col.nn <- c("none", "FIDIID", "PID_new", "MID_new", "Sex_new", "Phe_new")
col.t <- "ccccii"

ped <- ped %>%
  read_tsv(col_types = cols(`Family ID` = col_character(),
                            `Individual ID` = col_character(),
                            .default = col_skip())) %>%
  rename(FID_ped = `Family ID`, FID_new = `Individual ID`)

new_fam <- read_table2(new, col_names = col.nn, col_types = col.t) %>%
  separate(FIDIID, c("FID_new", "IID_new"), sep = "_",
           remove = F, extra = "merge") %>%
  left_join(ped, by = "FID_new") %>%
  mutate(IID_new_a = FID_new) %>%
  mutate(FID_new = ifelse(is.na(IID_new), FID_ped, FID_new)) %>%
  mutate(IID_new = ifelse(is.na(IID_new), IID_new_a, IID_new))

read_table2(old, col_names = col.n, col_types = col.t) %>%
  unite("FIDIID", FID, IID, sep = "_", remove = F) %>%
  right_join(new_fam, by = "FIDIID") %>%
  mutate(FID = ifelse(is.na(FID), FID_new, FID)) %>%
  mutate(IID = ifelse(is.na(IID), IID_new, IID)) %>%
  mutate(PID = ifelse(is.na(PID), PID_new, PID)) %>%
  mutate(MID = ifelse(is.na(MID), MID_new, MID)) %>%
  mutate(Sex = ifelse(is.na(Sex), Sex_new, Sex)) %>%
  mutate(Phe = ifelse(is.na(Phe), Phe_new, Phe)) %>%
  select(!!col.n) %>%
  write_tsv(out, col_names = F)

