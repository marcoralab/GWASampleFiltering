suppressMessages(library(readr))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
in_bim <- args[1]
out_bim <- args[2]
out_snplist <- args[3]

cat("in_bim: ", in_bim, "\n")
cat("out_bim: ", out_bim, "\n")
cat("out_snplist: ", out_snplist, "\n\n")

bim <- read_tsv(in_bim, col_names = F)
bim <- mutate(bim, X2 = paste(X1, X4, X6, X5, sep = ":"))

write_tsv(bim, out_bim, col_names = F)
write_tsv(select(bim, X2), out_snplist, col_names = F)
