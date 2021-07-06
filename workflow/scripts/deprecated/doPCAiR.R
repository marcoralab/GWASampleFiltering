#!/usr/bin/env Rscript

require(GENESIS)
require(GWASTools)
require(tibble)
require(readr)
require(dplyr)

fstem <- commandArgs(T)[1]
kingstem <- commandArgs(T)[2]

geno <- GenotypeData(GdsGenotypeReader(filename = paste0(fstem, ".gds")))

iids <- getScanID(geno)

KINGmat <- king2mat(file.kin0 = paste0(kingstem, ".kin0"),
                    file.kin = paste0(kingstem, ".kin"),
                    iids = iids)

FID <- getVariable(geno, "sample.annot/family")

PCA <- pcair(genoData = geno, kinMat = KINGmat, divMat = KINGmat)

eigenval <- PCA$values

namecols <- function(eigenvec) {
  colnames(eigenvec) <- paste0("PC", 1:ncol(eigenvec))
  return(eigenvec)
}

eigenvec <- PCA$vectors %>%
  namecols %>%
  as_tibble(rownames = "IID") %>%
  mutate(FID = FID) %>%
  select(FID, everything()) %>%
  select(FID, IID, paste0("PC", 1:10))

write_delim(eigenvec, paste0(fstem, ".eigenvec", delim = " "))
write_lines(eigenval, paste0(fstem, ".eigenval"))

