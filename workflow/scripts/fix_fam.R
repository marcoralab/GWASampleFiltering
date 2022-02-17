#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

if (!exists("snakemake")) {
  setClass("snakemake_fake", representation(
    params = "list", input = "list", output = "list",
    log = "list", wildcards = "list"))
  snakemake <- new("snakemake_fake",
    input = list(
      oldfam = "output/ADNI_3_chrall_pruned.fam",
      newfam = "output/ADNI_3_chrall_1kG_merged.fam",
      tgped = "reference/20130606_g1k.ped"
    ),
    params = list(),
    log = list(),
    output = list("output/ADNI_3_chrall_1kG_merged_fixed.fam"),
    wildcards = list()
  )
}

col.n <- c("FID", "IID", "PID", "MID", "Sex", "Phe")
col.nn <- c("none", "FIDIID", "PID_new", "MID_new", "Sex_new", "Phe_new")
col.t <- "ccccii"

new_fam_basic <- snakemake@input[["newfam"]] %>%
  read_table2(col_names = col.nn, col_types = col.t) %>%
  separate(FIDIID, c("FID_new", "IID_new"), sep = "_",
           remove = F, extra = "merge")

if (snakemake@input[["tgped"]]) {
  ped <- snakemake@input[["tgped"]] %>%
    read_tsv(col_types = cols(`Family ID` = col_character(),
                              `Individual ID` = col_character(),
                              .default = col_skip())) %>%
    rename(FID_ped = `Family ID`, FID_new = `Individual ID`)

    new_fam <- new_fam_basic %>%
      left_join(ped, by = "FID_new") %>%
      mutate(IID_new_a = FID_new) %>%
      mutate(FID_new = ifelse(is.na(IID_new), FID_ped, FID_new)) %>%
      mutate(IID_new = ifelse(is.na(IID_new), IID_new_a, IID_new))
} else {
  new_fam <- new_fam_basic
}

old_fam <- snakemake@input[["oldfam"]] %>%
  read_table2(col_names = col.n, col_types = col.t) %>%
  unite("FIDIID", FID, IID, sep = "_", remove = F)

left_join(new_fam, old_fam, by = "FIDIID") %>%
  mutate(FID = ifelse(is.na(FID), FID_new, FID)) %>%
  mutate(IID = ifelse(is.na(IID), IID_new, IID)) %>%
  mutate(PID = ifelse(is.na(PID), PID_new, PID)) %>%
  mutate(MID = ifelse(is.na(MID), MID_new, MID)) %>%
  mutate(Sex = ifelse(is.na(Sex), Sex_new, Sex)) %>%
  mutate(Phe = ifelse(is.na(Phe), Phe_new, Phe)) %>%
  select(!!col.n) %>%
  write_tsv(snakemake@output[["fixed"]], col_names = F)
