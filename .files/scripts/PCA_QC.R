#!/bin/bash

## ---- Load Required R packages ---- ##
library(tidyverse)  ## for data wrangling 

zscore = function(x){(x - mean(x)) / sd(x)}

## ---- File Input ---- ##
pca.file = commandArgs(TRUE)[1]
base_pops.file = commandArgs(TRUE)[2]
target_pops.file = commandArgs(TRUE)[3]
eigenval.file = commandArgs(TRUE)[4]
outfile = commandArgs(TRUE)[5]

##---- Read in Data ----##
# PCA file from plink
pca.raw <- as.tibble(read_delim(pca.file, delim = " ", col_names = F))

# population data file from 1000 genomes
base_pops.raw <- read_tsv(base_pops.file)

# population data from target set
target_pops.raw <- read_delim(target_pops.file, delim = " ", col_names = F)

# egien values 
eigenval.raw <- read_table(eigenval.file, col_names = F)

##---- Data wrangling ----##
##  rename column names in PCA file
colnames(pca.raw) <- c('FID', 'IID', paste0('PC', seq(1, ncol(pca.raw)-2, 1)))
##  standardize PC to a z-score
pca <- mutate_at(pca.raw, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10'), zscore)

##  Munge population dataframes from 1000 genomes
base_pops <- base_pops.raw %>% 
  dplyr:::rename(IID = `Individual ID`) %>%
  select(IID, Population) %>% 
  mutate(SuperPopulation = recode(base_pops.raw$Population, GBR = 'EUR', FIN = 'EUR', CHS = 'EAS', PUR = 'AMR', CDX = 'EAS', CLM = 'AMR', IBS = 'EUR', PEL = 'AMR', PJL = 'SAS', KHV = 'EAS', ACB = 'AFR', GWD = 'AFR', ESN = 'AFR', BEB = 'SAS', MSL = 'AFR', STU = 'SAS', ITU = 'SAS', CEU = 'EUR', YRI = 'AFR', CHB = 'EAS', JPT = 'EAS', LWK = 'AFR', ASW = 'AFR', MXL = 'AMR', TSI = 'EUR', GIH = 'SAS')) %>% 
  mutate(cohort = '1kgenomes') 

##  Munge target population dataframes
target_pops <- target_pops.raw %>% 
  select(X2) %>% 
  rename(IID = X2) %>% 
  mutate(Population = 'Sample', SuperPopulation = 'Sample', cohort = 'Sample', Population_2 = 'Sample') 

##  Munge PCA, base pop and target pop
pca <- target_pops %>% 
  bind_rows(base_pops) %>% 
  right_join(pca, by = 'IID')
 
##  Relevel Population Factor for ploting 
pca$Population_2 <- factor(pca$Population, level = c('Sample', 'GBR', 'FIN', 'IBS', 'CEU', 'TSI', 'PUR', 'CLM', 'PEL', 'MXL', 'CHS', 'CDX', 'KHV', 'CHB', 'JPT', 'PJL', 'BEB', 'STU', 'ITU', 'GIH', 'GWD', 'ESN', 'MSL', 'YRI', 'LWK', 'ASW', 'ACB')) 

## ---- Population Outliers ---- ##
##  For 1000 Genomes EUR population, calculate the mean and +- 6 standard dPCiations for each prinicipal compoent 
EurPca <- pca %>% 
  gather(key = 'PC', value = 'eigenvalue', PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
  filter(SuperPopulation == 'EUR') %>%
  group_by(SuperPopulation, PC) %>% 
  dplyr::summarize(mean = mean(eigenvalue), sd = sd(eigenvalue)) %>% 
  mutate(lower = mean - sd*6, upper = mean + sd*6) %>%
  mutate(PC = factor(PC, levels = paste0('PC', 1:10))) %>%
  arrange(PC)

EurPca <- as.data.frame(EurPca)

##  For eah invididule in the sample, determine weather samples are +- 6SD from the EUR population for each principal component
sample.pca <- pca %>% 
  filter(cohort == 'Sample') %>%
  mutate(PC1.outliers = PC1 < EurPca[EurPca$PC == 'PC1', 'lower'] | PC1 > EurPca[EurPca$PC == 'PC1', 'upper']) %>% 
  mutate(PC2.outliers = PC2 < EurPca[EurPca$PC == 'PC2', 'lower'] | PC2 > EurPca[EurPca$PC == 'PC2', 'upper']) %>% 
  mutate(PC3.outliers = PC3 < EurPca[EurPca$PC == 'PC3', 'lower'] | PC3 > EurPca[EurPca$PC == 'PC3', 'upper']) %>% 
  mutate(PC4.outliers = PC4 < EurPca[EurPca$PC == 'PC4', 'lower'] | PC4 > EurPca[EurPca$PC == 'PC4', 'upper']) %>% 
  mutate(PC5.outliers = PC5 < EurPca[EurPca$PC == 'PC5', 'lower'] | PC5 > EurPca[EurPca$PC == 'PC5', 'upper']) %>% 
  mutate(PC6.outliers = PC6 < EurPca[EurPca$PC == 'PC6', 'lower'] | PC6 > EurPca[EurPca$PC == 'PC6', 'upper']) %>%
  mutate(PC7.outliers = PC7 < EurPca[EurPca$PC == 'PC7', 'lower'] | PC7 > EurPca[EurPca$PC == 'PC7', 'upper']) %>%
  mutate(PC8.outliers = PC8 < EurPca[EurPca$PC == 'PC8', 'lower'] | PC8 > EurPca[EurPca$PC == 'PC8', 'upper']) %>%
  mutate(PC9.outliers = PC9 < EurPca[EurPca$PC == 'PC9', 'lower'] | PC9 > EurPca[EurPca$PC == 'PC9', 'upper']) %>%
  mutate(PC10.outliers = PC10 < EurPca[EurPca$PC == 'PC10', 'lower'] | PC10 > EurPca[EurPca$PC == 'PC10', 'upper']) %>%
  mutate(pop.outliers =  PC1.outliers == T | PC2.outliers == TRUE | PC3.outliers == TRUE | PC4.outliers == TRUE | PC5.outliers == TRUE | PC6.outliers == TRUE | PC7.outliers == TRUE | PC8.outliers == TRUE | PC9.outliers == TRUE  | PC10.outliers == TRUE)

##	---- Write out outliers ---- ##
exclude_pop_outliers <- sample.pca %>% 
  filter(pop.outliers == TRUE) %>% 
  select(FID, IID)

write_tsv(exclude_pop_outliers, outfile, col_names = F)




































