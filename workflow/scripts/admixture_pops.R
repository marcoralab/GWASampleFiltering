library(readr)
library(dplyr)

# .Pop file for admixture
pop_info <- snakemake@input[["spop"]] %>%
  read_tsv(col_types = cols(.default = "c"))
  
refpops <- snakemake@input[['pops']] %>%
  read_table2(col_types = "ccc") %>%
  rename(pop = Population) %>%
  left_join(pop_info, by = "pop")

message('\nWriting .pop file: ', snakemake@output[[1]])
snakemake@input[["fam"]] %>%
  read_tsv(col_names = c("FID", "IID"), col_types = "cc----") %>%
  left_join(refpops, by = c("FID", "IID")) %>%
  mutate(spop = ifelse(is.na(spop), "-", spop)) %>%
  pull(spop) %>%
  write_lines(snakemake@output[[1]])
