#!/bin/bash
## Heterozygosity check for GWAS data

## ---- Load Pacakges ---- ##
library(tidyverse)
library(ggplot2)

het.file = commandArgs(TRUE)[1]
outfile = commandArgs(TRUE)[2]

het.raw <- as.tibble(read.table(het.file, header = TRUE, check.names = FALSE, as.is = TRUE, colClasses = c("character","character","numeric","numeric","numeric","numeric")))

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



