#!/bin/bash
##	R script for identifing samples that are crypticaly related from PLINK files.

require(tidyverse, quietly = T)

##  Arguments: 
# 1: .genome output file from PLINK
# 2: .fam file from PLINK containing subject information
# 3: .txt file to be outputed
genome.file = commandArgs(TRUE)[1]
fam.file = commandArgs(TRUE)[2]
outfile = commandArgs(TRUE)[3]

##  read in .genome file
dat.inter <- as.tibble(read.table(genome.file, header = TRUE, check.names = FALSE, as.is = TRUE))

##  select samples with kinship cofficents > 0.1
ibdcoeff <- filter(dat.inter, PI_HAT > 0.1)

##  Iterativly remove subjects with the highest number of pairwise kinship cofficents > 0.1 
##  see http://www.stat-gen.org/tut/tut_preproc.html 
related.samples <- NULL
while(nrow(ibdcoeff) >0 ){
  sample.counts <- arrange(plyr:::count(c(ibdcoeff$IID1, ibdcoeff$IID2)), -freq)
  rm.sample <- sample.counts[1,'x']
  cat("Removing sample", as.character(rm.sample), 'to closely related to', sample.counts[1, 'freq'], 'other samples. \n')
  
  ibdcoeff <- ibdcoeff[ibdcoeff$IID1 != rm.sample & ibdcoeff$IID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

##  Read in PLINK .fam file, extract FID and IID of samples to be excluded
fam <- read_delim(fam.file, delim = " ", col_names = F)
exclude.samples <- fam %>% 
  filter(X2 %in% related.samples) %>% 
  select(X1, X2) %>% 
  rename(FID = X1, IID = X2)

##  write out samples to be excluded
write_tsv(exclude.samples, outfile, col_names = T)


