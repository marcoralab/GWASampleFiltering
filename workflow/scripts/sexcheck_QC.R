#!/usr/bin/env Rscript
## File for identifing sex discordance

require(tidyverse, quietly = T)

##  Arguments:
# .sexcheck output file from PLINK
sexcheck_file <- snakemake@input[[1]]
# .txt file to be outputed
outfile <- snakemake@output[[1]]
if (is.null(outfile)) {print("Outfile is required! \n")}

sexcheck <- sexcheck_file %>%
  read.table(header = T, check.names = F, as.is = T,
             colClasses = c(
               "character", "character", "integer",
               "integer", "character", "numeric")) %>%
  as_tibble()

exclude_samples <- sexcheck %>%
  filter(STATUS == "PROBLEM") %>%
  select(FID, IID)

cat("Removing sample",
    as.character(exclude_samples$IID),
    "due to sex discordance. \n")

##  write out samples to be excluded
write_tsv(exclude_samples, outfile, col_names = T)
