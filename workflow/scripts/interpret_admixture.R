#!/usr/bin/env Rscript
## ========================================================================== ##
## Assign ancestry for admixture output
## ========================================================================== ##

suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)

if (!exists("snakemake")) {
  cohort_name <- "ADSP17kWGS"
  inputs <- lapply(
    list(Qraw = "output/%s_1kG_filtered.3.Q",
         fam = "output/%s_1kG_admixfiltered.fam",
         pops = "output/%s_1kG_admixfiltered.pop",
         pcs = "output/%s_cluster_pops.tsv"),
    sprintf,
    cohort_name)
  output = sprintf("output/%s_1kG_merged_interpret_admixture.tsv", cohort_name)
  setClass("snakemake_fake", representation(
    params = "list", input = "list", output = "list",
    log = "list", wildcards = "list"))
  snakemake <- new("snakemake_fake",
    input = inputs,
    params = list(),
    log = list(),
    output = list(output),
    wildcards = list()
  )
}

## Infiles
Qraw = snakemake@input[["Qraw"]]
path_famfile = snakemake@input[["fam"]]
path_popfile = snakemake@input[["pops"]]
path_pcs = snakemake@input[["pcs"]]

## Outfile
supervised_assign = snakemake@output[[1]]

# Fam and popfiles
## ======================================##
message("Reading pop file \n")
popfile <- path_popfile %>%
  read_table(col_names = c("super_pop"), col_types = "c")

super_labels_frompop <- popfile %>%
  distinct %>%
  filter(super_pop != "-") %>%
  mutate(col = paste0("k", 1:nrow(.))) %>%
  deframe

message("Reading fam fixed file \n")
famfile <- path_famfile %>%
  read_table2(col_names = c("FID", "IID"), col_types = "cc----") %>%
  bind_cols(popfile)

message("Reading pcs file \n")
pcs <- path_pcs %>%
  read_tsv(col_types = cols(
    .default = "d",
    superpop_infered = "c",
    FID = "c",
    IID = "c",
    Population = "c",
    superpop = "c",
    cohort = "c")) %>%
  rename(pca_super_pop = superpop_infered)

# Interpreting supervised admixture output #
## ======================================##
message("Reading supervised admixture output \n")
tbl_super <- Qraw %>%
  read_table(col_names = F, col_types = cols(.default = "d")) %>%
  bind_cols(famfile) %>%
  left_join(select(pcs, FID, IID, pca_super_pop), by = c('FID', 'IID')) %>%
  filter(super_pop == "-") %>%
  select(-super_pop) %>%
  rename_with( ~ str_replace(.x, "^X", "k"))

super_labels <- tbl_super %>%
  group_by(pca_super_pop) %>%
  summarise(across(starts_with("k"), mean), .groups = "drop") %>%
  mutate(max = pmap_chr(select(., starts_with("k")),
                        ~ c(...) %>% which.max %>% names)) %>%
  arrange(max) %>%
  select(pca_super_pop, max) %>%
  deframe

if (length(super_labels) < length(super_labels_frompop)) {
  test_max <- enframe(super_labels) %>% rename(max = value)
  test_pop <- enframe(super_labels_frompop) %>% rename(pop = value)
  test_names <- full_join(test_max, test_pop, by = "name")
  if (!all(test_names$pop == test_names$max, na.rm=T)) {
    print(test_names)
    stop("Unsure of population names!")
  }
  super_labels <- super_labels_frompop
}

out_start <- tbl_super %>%
  rename(!!!super_labels) %>%
  relocate(FID, IID)

super_guess <- out_start %>%
  bind_rows(tibble(
    EUR = 0, EAS = 0, SAS = 0, AFR = 0, AMR = 0, dummy = T)) %>%
  mutate(
    across(where(is.double), ~ replace_na(.x, 0)),
    admixture_super_pop_max = pmap_chr(
      select(., !!names(super_labels)),
      ~ c(...) %>% which.max %>% names),
    admixture_super_pop = case_when(
      (EUR > 0.85 & EAS < 0.1 & SAS < 0.1 & AFR < 0.1 & AMR < 0.1) ~ "EUR",
      (EAS > 0.51) ~ "EAS",
      (SAS > 0.51) ~ "SAS",
      (AFR > 0.3 & EAS < 0.1 & SAS < 0.1 & AFR > AMR) ~ "AFR",
      (AMR > 0.1 & EAS < 0.1 & SAS < 0.1 ) ~ "AMR",
      TRUE ~ "Other")) %>%
  filter(is.na(dummy)) %>%
  select(FID, IID, admixture_super_pop_max, admixture_super_pop)

out_super <- left_join(out_start, super_guess, by = c("FID", "IID"))

message("Exporting out to ", supervised_assign, "\n")
out_super %>% write_tsv(supervised_assign)
