#!/usr/bin/env Rscript
##	R script for concating sample exlusion lists from PLINK GWAS QC.

require(dplyr, quietly = T, warn.conflicts = F)
require(readr, quietly = T, warn.conflicts = F)


# 1: call rate fail sample list
irem.file <- snakemake@input[['SampCallRate']]
# 2: heterozygocity fail sample list
het.file <- snakemake@input[['het']]
# 3: discordent sex sample list
sex.file <- snakemake@input[['sex']]
# 4: population outliers sample list
pca.file <- snakemake@input[['pca']]
# 5: duplicates/relatives sample list
rel.file <- snakemake@input[['relat']]
# 6: .txt file to be outputed
outfile <- snakemake@output[['out']]
outfile_d <- snakemake@output[['out_distinct']]

##  ---- Read in sample exclusion files files ---- ##
irem.raw <- read_tsv(irem.file, col_names = F, col_types = 'cc') %>%
  mutate(why = 'mis')
het.raw <- read_tsv(het.file, col_names = T, col_types = 'cc') %>%
  mutate(why = 'het')
sex.raw <- read_tsv(sex.file, col_names = T, col_types = 'cc') %>%
  mutate(why = 'sex')
pca.raw <- read_tsv(pca.file, col_names = F, col_types = 'cc') %>%
  mutate(why = 'pca')
rel.raw <- read_tsv(rel.file, col_names = T, col_types = 'cc') %>%
  mutate(why = 'rel')

##  ---- Data wrangling ---- ##
## IF .irem file is empty, make empty tibble
fill_blanks <- function (df) {
  if (nrow(df) == 0) {
    out <- tibble(FID = NA, IID = NA)
  } else if ('X1' %in% names(df)) {
    out <- df %>%
      rename(FID = X1, IID = X2)
  } else {
    out <- df
  }
  return(out)
}

irem <- fill_blanks(irem.raw)
het <- fill_blanks(het.raw)
sex <- fill_blanks(sex.raw)
pca <- fill_blanks(pca.raw)
rel <- fill_blanks(rel.raw)

excluded <- het %>%
  bind_rows(irem) %>%
  bind_rows(sex) %>%
  bind_rows(pca) %>%
  bind_rows(rel) %>%
  filter(!is.na(FID))
message("Excluded:")
print(data.frame(excluded))

##  write out samples to be excluded
write_tsv(excluded, outfile, col_names = T)

excluded %>%
  distinct %>%
  write_tsv(outfile_d, col_names = T)
