suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)
InBim <- args[1]
OutBim <- args[2]
OutSnplist <- args[3]

cat('InBim: ', InBim, '\n')
cat('OutBim: ', OutBim, '\n')
cat('OutSnplist: ', OutSnplist, '\n\n')

bim <- read_tsv(InBim, col_names = F) 
bim <- mutate(bim, X2 = paste(X1, X4, X5, X6, sep = ":"))

write_tsv(bim, OutBim, col_names = F)
write_tsv(select(bim, X2), OutSnplist, col_names = F)
