suppressMessages(library(readr))
suppressMessages(library(dplyr))

in_bim = snakemake@input[['flipped']]
out_bim = snakemake@output[['bim']]
out_snplist = snakemake@output[['snplist']]

cat("in_bim: ", in_bim, "\n")
cat("out_bim: ", out_bim, "\n")
cat("out_snplist: ", out_snplist, "\n\n")

bim <- read_tsv(in_bim, col_names = F,
  col_types = cols(
  X1 = col_double(),
  X2 = col_character(),
  X3 = col_double(),
  X4 = col_double(),
  X5 = col_character(),
  X6 = col_character()
))
bim <- mutate(bim, X2 = paste(X1, X4, X6, X5, sep = ":"))

write_tsv(bim, out_bim, col_names = F)
write_tsv(select(bim, X2), out_snplist, col_names = F)
