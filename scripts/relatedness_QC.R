#!/bin/bash
##	R script for identifing samples that are crypticaly related from PLINK files.

suppressMessages(require(tidyverse))
suppressMessages(require(Hmisc))
suppressMessages(require(magrittr))


##  Arguments:
# 1: .genome output file from PLINK
# 2: pi_hat threshold (0.1875 recommended)
# 3: Family based study, T/F
# 4: .tsv file to be outputed
genome.file <- commandArgs(TRUE)[1]
threshold <- as.numeric(commandArgs(TRUE)[2])
Family <- as.logical(commandArgs(TRUE)[3])
king <- as.logical(commandArgs(TRUE)[4])
outfile <- commandArgs(TRUE)[5]
rdat <- commandArgs(TRUE)[6]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Read in Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

if (king) {
  dat.inter <- read_table2(genome.file, col_types = cols(
    .default = col_double(),
    FID1 = col_character(),
    ID1 = col_character(),
    FID2 = col_character(),
    ID2 = col_character(),
    N_SNP = col_integer(),
    InfType = col_character()
  )) %>%
    rename(IID1 = ID1, IID2 = ID2, PI_HAT = PropIBD)
} else {
  dat.inter <- read_table2(genome.file, col_types = cols(
    .default = col_double(),
    FID1 = col_character(),
    IID1 = col_character(),
    FID2 = col_character(),
    IID2 = col_character(),
    RT = col_character()
  ))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  IBD relationship
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# IBD relationship table
# https://github.com/WheelerLab/GWAS_QC/blob/master/example_pipelines/QC%20Analysis%20-%20Cox%20Lab%20Projects.pdf
rel_tab <- tibble(relationship = c("unrelated", "identical-twins",
                                   "parent-child", "full-siblings",
                                   "half-siblings", "grandparent-grandchild",
                                   "avuncular", "half-avuncular",
                                   "first-cousin", "half-first-cousin",
                                   "half-sibling-first-cousin"),
  pi_hat = c(0, 1, 0.5, 0.5, 0.25, 0.25, 0.25, 0.125, 0.125, 0.0625, 0.375),
  z0 = c(1, 0, 0, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.875, 0.375),
  z1 = c(0, 0, 1, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.125, 0.5),
  z2 = c(0, 1, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.125)
)
##  Remove duplicate values, assign duplicate lables to same lable
dup_relationships <- c("grandparent-grandchild", "avuncular", "half-avuncular")
rel_tab_filt <- rel_tab %>%
  filter(relationship %nin% dup_relationships) %>%
  mutate(relationship = ifelse(relationship == "half-siblings", "2nd degree",
                               ifelse(relationship == "first-cousin", "3rd degree", relationship)))

# Match Pi-Hat, Z0, Z1, Z2 to closet values in IBD relationship table

closest <- function(vals, ref) {
  fc <- Vectorize(function (x) {
    ref[which.min(abs(ref - x))]
  }) #finds closest
  fc(vals)
}



if (!king) {
  dat.inter %<>%
    mutate(pi_hat = closest(PI_HAT, rel_tab_filt$pi_hat),
           z0 = closest(Z0, rel_tab_filt$z0),
           z1 = closest(Z1, rel_tab_filt$z1),
           z2 = closest(Z2, rel_tab_filt$z2)) %>%
    left_join(rel_tab_filt, by = c("pi_hat", "z0", "z1", "z2")) %>%
    mutate(relationship = ifelse(is.na(relationship), "OA", relationship))
} else {
  dat.inter %<>%
    mutate(relationship = InfType) %>%
    mutate(relationship = ifelse(relationship %in% c("1st", "2nd", "3rd", "4th"),
                                 paste(relationship, "degree"),
                                 relationship)) %>%
    mutate(relationship = ifelse(relationship == "UN",
                                 "unrelated",
                                 relationship)) %>%
    mutate(relationship = ifelse(relationship == "Dup/MZ",
                                 "identical-twins",
                                 relationship)) %>%
    mutate(relationship = ifelse(relationship == "PO",
                                 "parent-child",
                                 relationship)) %>%
    mutate(relationship = ifelse(relationship == "FS",
                                 "full-siblings",
                                 relationship))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Exclude Samples
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# select samples with kinship cofficents > 0.1875
# https://link.springer.com/protocol/10.1007/978-1-60327-367-1_19

if (Family == F | Family == "F"){
  ibdcoeff <- dat.inter %>%
    filter(PI_HAT > threshold)
} else {
  ibdcoeff <- dat.inter %>%
    filter(relationship == "identicial-twins")
}

if (king) {
  ibd_tab <- ibdcoeff %>%
    select(FID1, IID1, FID2, IID2, IBS0, Kinship, PI_HAT, relationship)
} else {
  ibd_tab <- ibdcoeff %>%
    select(FID1, IID1, FID2, IID2, Z0, Z1, Z2, PI_HAT, relationship)
}

# Iterativly remove subjects with the highest number of pairwise kinship cofficents > threshold
# see http://www.stat-gen.org/tut/tut_preproc.html

if (Family == F | Family == "F") {
  message("Working with unrelated samples.")
  ibdcoeff %<>%
    mutate(FI1 = paste0(FID1, "🤣", IID1), FI2 = paste0(FID2, "🤣", IID2))
  related.samples <- NULL
  excluded <- c()
  fam_table <- tibble(FID = c("deleteme"), IID = c("deleteme"), Related = c("deleteme"))
  while (nrow(ibdcoeff) > 0 ) {
    sample.counts <- arrange(plyr:::count(c(ibdcoeff$FI1, ibdcoeff$FI2)), -freq)
    rm.sample <- sample.counts[1, "x"]
    ID <- str_split(rm.sample, "🤣")[[1]]
    FID <- ID[1]
    IID <- ID[2]
    remtxt <- sprintf("closely related to %i other samples.",
                      sample.counts[1, "freq"])
    message(paste("Removing sample", IID, remtxt))
    ft <- tibble(FID = FID, IID = IID, Related = remtxt)
    fam_table <- fam_table %>%
      bind_rows(ft)
    ibdcoeff <- ibdcoeff[ibdcoeff$FI1 != rm.sample &
                         ibdcoeff$FI2 != rm.sample, ]
    related.samples <- c(as.character(rm.sample), related.samples)
  }
  # Don't simplify fam_table. It will break.
  fam_table <- fam_table %>%
    filter(Related != "deleteme")
  exclude.samples <- tibble(FI = as.character(related.samples)) %>%
    separate(FI, c("FID", "IID"), sep = "🤣")
} else {
  message("Working with related samples.")
  fam_table <- as.data.frame(ibdcoeff)
  related.samples <- c(ibdcoeff$IID1, ibdcoeff$IID2)
  exclude.samples <- tibble(FID = character(), IID = character())
}

##  write out samples to be excluded
write_tsv(exclude.samples, outfile, col_names = T)
save(dat.inter, rel_tab, fam_table, ibd_tab, king, file = rdat)
