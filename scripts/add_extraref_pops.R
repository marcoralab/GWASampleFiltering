#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)

args <- commandArgs(T)

refpops <- args[1]
extravcf <- args[2]
subpop <- args[3]
out <- args[4]
out2 <- args[5]

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
