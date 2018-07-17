#!/bin/bash
##	R script for identifing samples that are crypticaly related from PLINK files.

suppressMessages(require(tidyverse))
suppressMessages(require(Hmisc))


##  Arguments: 
# 1: .genome output file from PLINK
# 2: .fam file from PLINK containing subject information
# 3: Family based study, T/F
# 4: .txt file to be outputed
genome.file = commandArgs(TRUE)[1]
fam.file = commandArgs(TRUE)[2]
Family = commandArgs(TRUE)[3]
outfile = commandArgs(TRUE)[4]

##----------------------------##
##  Read in Data
##----------------------------##

dat.inter.raw <- as.tibble(read.table(genome.file, header = TRUE, check.names = FALSE, as.is = TRUE))
dat.inter <- dat.inter.raw

##----------------------------##
##  IBD relationship
##----------------------------##

# IBD relationship table
# https://github.com/WheelerLab/GWAS_QC/blob/master/example_pipelines/QC%20Analysis%20-%20Cox%20Lab%20Projects.pdf
rel_tab <- tibble(relationship = c('unrelated', 'identicial-twins', 'parent-child', 'full-siblings', 'half-siblings', 'grandparent-grandchild', 'avuncular', 'half-avuncular', 'first-cousin', 'half-first-cousin', 'half-sibling-first-cousin'),
                  pi_hat = c(0, 1, 0.5, 0.5, 0.25, 0.25, 0.25, 0.125, 0.125, 0.0625, 0.375),
                  z0 = c(1, 0, 0, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.875, 0.375),
                  z1 = c(0, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.125, 0.5),
                  z2 = c(0, 1, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.125)
)
##  Remove duplicate values, assign duplicate lables to same lable
rel_tab_filt <- filter(rel_tab, relationship %nin% c('grandparent-grandchild', 'avuncular', 'half-avuncular'))
rel_tab_filt[rel_tab_filt$relationship == 'half-siblings', ]$relationship <- "half-siblings \n grandparent-grandchild \n avuncular"
rel_tab_filt[rel_tab_filt$relationship == 'first-cousin', ]$relationship <- "first-cousin \n half-avuncular"

## Match Pi-Hat, Z0, Z1, Z2 to closet values in IBD relationship table
dat.inter$pi_hat <- sapply(dat.inter$PI_HAT, function(x){
  out <- rel_tab$pi_hat[which(abs(rel_tab$pi_hat-x)==min(abs(rel_tab$pi_hat-x)))]
  out[1]
})
dat.inter$z0 <- sapply(dat.inter$Z0, function(x){
  out <- rel_tab$z0[which(abs(rel_tab$z0-x)==min(abs(rel_tab$z0-x)))]
  out[1]
})
dat.inter$z1 <- sapply(dat.inter$Z1, function(x){
  out <- with(rel_tab, z1[which(abs(z1-x)==min(abs(z1-x)))])
  out[1]
})
dat.inter$z2 <- sapply(dat.inter$Z2, function(x){
  out <- rel_tab$z2[which(abs(rel_tab$z2-x)==min(abs(rel_tab$z2-x)))]
  out[1]
})

## Left join 
dat.inter <- dat.inter %>% left_join(rel_tab_filt, by = c('pi_hat', 'z0', 'z1', 'z2')) 
dat.inter <- mutate(dat.inter, relationship = ifelse(is.na(relationship), 'OA', relationship))

##----------------------------##
##  Exclude Samples
##----------------------------##

##  select samples with kinship cofficents > 0.1
if(Family == F){
  ibdcoeff <- dat.inter %>% 
    filter(PI_HAT > 0.1) %>% 
    select(FID1, IID1, FID2, IID2, Z0, Z1, Z2, PI_HAT, relationship)
} else {
  ibdcoeff <- dat.inter %>% 
    filter(relationship == 'identicial-twins') %>% 
    select(FID1, IID1, FID2, IID2, Z0, Z1, Z2, PI_HAT, relationship)
}


##  Iterativly remove subjects with the highest number of pairwise kinship cofficents > 0.1 
##  see http://www.stat-gen.org/tut/tut_preproc.html 
if(Family == F){
  related.samples <- NULL
  while(nrow(ibdcoeff) >0 ){
    sample.counts <- arrange(plyr:::count(c(ibdcoeff$IID1, ibdcoeff$IID2)), -freq)
    rm.sample <- sample.counts[1,'x']
    cat("Removing sample", as.character(rm.sample), 'to closely related to', sample.counts[1, 'freq'], 'other samples. \n')
    
    ibdcoeff <- ibdcoeff[ibdcoeff$IID1 != rm.sample & ibdcoeff$IID2 != rm.sample,]
    related.samples <- c(as.character(rm.sample), related.samples)
  }
} else {
  related.samples <- c(ibdcoeff$IID1, ibdcoeff$IID2)
}

##  Read in PLINK .fam file, extract FID and IID of samples to be excluded
fam <- read_delim(fam.file, delim = " ", col_names = F)
exclude.samples <- fam %>% 
  filter(X2 %in% related.samples) %>% 
  select(X1, X2) %>% 
  rename(FID = X1, IID = X2)

##  write out samples to be excluded
write_tsv(exclude.samples, outfile, col_names = T)


