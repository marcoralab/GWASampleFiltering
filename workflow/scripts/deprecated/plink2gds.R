#!/usr/bin/env Rscript

require(SNPRelate)

pstem <- commandArgs(T)[1]

snpgdsBED2GDS(bed.fn = paste0(pstem, ".bed"),
              bim.fn = paste0(pstem, ".bim"),
              fam.fn = paste0(pstem, ".fam"),
              out.gdsfn = paste0(pstem, ".gds"),
              family = T)

