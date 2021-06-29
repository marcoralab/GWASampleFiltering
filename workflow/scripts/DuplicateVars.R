suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))

Infile <- snakemake@input[[1]]
Outfile <- snakemake@output[[1]]

dupvar <- read_tsv(Infile, col_names = T)
dupvar %>%
  separate(IDS, c('Var1', 'Var2', 'Var3'), sep = " ") %>%
  select(Var1) %>%
  write_tsv(Outfile)
