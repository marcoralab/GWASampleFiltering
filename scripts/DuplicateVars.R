suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))

args = commandArgs(trailingOnly=TRUE)
Infile <- args[1]

dupvar <- read_tsv(Infile, col_names = T)
dupvar %>%
  separate(IDS, c('Var1', 'Var2', 'Var3'), sep = " ") %>%
  select(Var1) %>%
  write_tsv(paste0(Infile, '.delete'))
