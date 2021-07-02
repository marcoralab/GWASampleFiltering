#!/bin/bash
##	File for identifing sex discrodance
require(tidyverse, quietly =T)

##  Arguments:
# 1: .sexcheck output file from PLINK
# 3: .txt file to be outputed
sexcheck.file = snakemake@input[[1]]
outfile = snakemake@output[[1]]
if(is.null(outfile)){print('Outfile is required! \n')}

sexcheck <- as.tibble(read.table(sexcheck.file, header = TRUE, check.names = FALSE, as.is = TRUE,
                                 colClasses = c("character","character","integer","integer","character","numeric")))
exclude.samples <- sexcheck %>%
  filter(STATUS == 'PROBLEM') %>%
  select(FID, IID)

cat("Removing sample", as.character(exclude.samples$IID), 'due to sex discordance. \n')

##  write out samples to be excluded
write_tsv(exclude.samples, outfile, col_names = T)