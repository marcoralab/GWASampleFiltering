#!/bin/bash
##	R script for concating sample exlusion lists from PLINK GWAS QC.

require(tidyverse, quietly = T)

##  Arguments: 
# 1: call rate fail sample list
# 2: heterozygocity fail sample list
# 3: discordent sex sample list
# 4: population outliers sample list
# 5: duplicates/relatives sample list
# 6: .txt file to be outputed

irem.file = commandArgs(TRUE)[1]
het.file = commandArgs(TRUE)[2]
sex.file = commandArgs(TRUE)[3]
pca.file = commandArgs(TRUE)[4]
rel.file = commandArgs(TRUE)[5]
outfile = commandArgs(TRUE)[6]

##  ---- Read in sample exclusion files files ---- ##
irem.raw <- read_tsv(irem.file, col_names = F)
het.raw <- read_tsv(het.file, col_names = T)
sex.raw <- read_tsv(sex.file, col_names = T)
pca.raw <- read_tsv(pca.file, col_names = F)
rel.raw <- read_tsv(rel.file, col_names = T)

##  ---- Data wrangling ---- ##

irem <- irem.raw %>% 
  rename(FID = X1, IID = X2) 

het <- het.raw 

sex <- sex.raw 

pca <- pca.raw %>% 
  rename(FID = X1, IID = X2) 
  
rel <- rel.raw 

excluded <- het %>% 
  full_join(irem, by = c('FID', 'IID')) %>%
  full_join(sex, by = c('FID', 'IID')) %>% 
  full_join(pca, by = c('FID', 'IID')) %>%
  full_join(rel, by = c('FID', 'IID'))

##  write out samples to be excluded
write_tsv(excluded, outfile, col_names = T)
