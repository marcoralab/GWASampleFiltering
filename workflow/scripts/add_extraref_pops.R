#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)

refpops = snakemake@input[['main']]
extravcf = snakemake@input[['extra']]
subpop = snakemake@params[['extra_ref_code']]
out = snakemake@output[[1]]
out2 = snakemake@output[[2]]

extraref <- sprintf('zgrep -m1 "^#CHROM" %s', extravcf) %>%
  pipe() %>%
  readLines() %>%
  strsplit("\t") %>%
  {tibble(FIDIID = .[[1]][- (1:9)])} %>%
  separate(FIDIID, c("FID", "IID"), sep = "_", extra = "merge") %>%
  mutate(Population = subpop) %>%
  bind_rows(read_delim(refpops, delim = " ", col_types = "ccc"), .)

write_delim(extraref, out, delim = " ")

extraref %>% distinct(Population) %>% pull(Population) %>% write_lines(out2)
