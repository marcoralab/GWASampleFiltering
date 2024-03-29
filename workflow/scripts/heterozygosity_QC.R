#!/usr/bin/env Rscript

## Heterozygosity check for GWAS data

## ---- Load Pacakges ---- ##
library(dplyr)
library(tibble)
library(readr)

het.file = snakemake@input[[1]]
outfile = snakemake@output[[1]]

het.raw <- as_tibble(read.table(het.file, header = TRUE, check.names = FALSE, as.is = TRUE, colClasses = c("character","character","numeric","numeric","numeric","numeric")))

## caluclate heterozygosity
het <- het.raw %>%
  rename(O = `O(HOM)`, E = `E(HOM)`, N = `N(NM)`) %>%
  mutate(Het = (N - O) / N)

##  Calculate exclusion thresholds
upper.het <- mean(het$Het) + sd(het$Het)*3
lower.het <- mean(het$Het) - sd(het$Het)*3

##  Exclusion of samples
het <- het %>%
  mutate(exclude = ifelse(Het >= upper.het | Het <= lower.het, TRUE, FALSE))

exclude.samples <- het %>% filter(exclude == TRUE)

sample.out <- exclude.samples %>%
  select(FID, IID)

cat("Removing sample", as.character(exclude.samples$IID), 'due to outlying heterozygosity. \n')

##  write out samples to be excluded
write_tsv(sample.out, outfile, col_names = T)
