#!/usr/bin/env Rscript
##	R script for identifing samples that are crypticaly related from PLINK files.

suppressMessages(require(readr))
suppressMessages(require(dplyr))
suppressMessages(require(tibble))
suppressMessages(require(stringr))
suppressMessages(require(tidyr))
suppressMessages(require(magrittr))

`%nin%` <- Negate(`%in%`)

genome.file <- snakemake@params[['geno']]
threshold <- as.numeric(snakemake@params[['threshold']])
Family <- as.logical(snakemake@params[['Family']])
king <- as.logical(snakemake@params[['king']])
outfile <- snakemake@output[['out']]
rdat <- snakemake@output[['rdat']]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Read in Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

if (king) {
  kin0 <- F
  kin <- F
  dat.inter.kin0 <- paste0(genome.file, ".kin0")
  dat.inter.kin <- paste0(genome.file, ".kin")
  if (file.exists(dat.inter.kin0)) {
    kin0 <- T
    dat.inter.kin0 %<>% read_table2(col_types = cols(
      .default = col_double(),
      FID1 = col_character(),
      ID1 = col_character(),
      FID2 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer(),
      InfType = col_character()
    )) %>%
      rename(IID1 = ID1, IID2 = ID2, PI_HAT = PropIBD)
    dat.inter.all.kin0 <- read_table2(paste0(genome.file, ".all.kin0"), col_types = cols(
      .default = col_double(),
      FID1 = col_character(),
      ID1 = col_character(),
      FID2 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
    )) %>%
      rename(IID1 = ID1, IID2 = ID2) %>%
      mutate(PI_HAT = ifelse(Kinship > 0, 2*Kinship, 0))
  }
  if (file.exists(dat.inter.kin)) {
    kin <- T
    dat.inter.kin <- read_table2(paste0(genome.file, ".kin"), col_types = cols(
      .default = col_double(),
      FID = col_character(),
      ID1 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer(),
      InfType = col_character()
    )) %>%
      rename(IID1 = ID1, IID2 = ID2, FID1 = FID, PI_HAT = PropIBD) %>%
      mutate(FID2 = FID1)
    dat.inter.all.kin <- read_table2(paste0(genome.file, ".all.kin"), col_types = cols(
      .default = col_double(),
      FID = col_character(),
      ID1 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
    )) %>%
      rename(IID1 = ID1, IID2 = ID2, FID1 = FID) %>%
      mutate(FID2 = FID1, PI_HAT = ifelse(Kinship > 0, 2*Kinship, 0))
  }
  if (kin0 & kin) {
    dat.inter.all <- bind_rows(dat.inter.all.kin0, dat.inter.all.kin) %>%
      distinct(FID1, IID1, FID2, IID2, .keep_all = T)
    dat.inter <- bind_rows(dat.inter.kin0, dat.inter.kin) %>%
      distinct(FID1, IID1, FID2, IID2, .keep_all = T)
  } else if (kin0) {
    dat.inter.all <- dat.inter.all.kin0
    dat.inter <- dat.inter.kin0
  } else if (kin) {
    dat.inter.all <- dat.inter.all.kin
    dat.inter <- dat.inter.kin
  }
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
    mutate(FI1 = paste0(FID1, "_-_-tempsep-_-_", IID1), FI2 = paste0(FID2, "_-_-tempsep-_-_", IID2))
  related.samples <- NULL
  excluded <- c()
  fam_table <- tibble(FID = c("deleteme"), IID = c("deleteme"), Related = c("deleteme"))
  while (nrow(ibdcoeff) > 0 ) {
    sample.counts <- arrange(plyr:::count(c(ibdcoeff$FI1, ibdcoeff$FI2)), -freq)
    rm.sample <- sample.counts[1, "x"]
    ID <- str_split(rm.sample, "_-_-tempsep-_-_")[[1]]
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
    separate(FI, c("FID", "IID"), sep = "_-_-tempsep-_-_")
} else {
  message("Working with related samples.")
  fam_table <- as.data.frame(ibdcoeff)
  related.samples <- c(ibdcoeff$IID1, ibdcoeff$IID2)
  exclude.samples <- tibble(FID = character(), IID = character())
}

##  write out samples to be excluded
write_tsv(exclude.samples, outfile, col_names = T)
if (king) {
  save(dat.inter, dat.inter.all, rel_tab, fam_table, ibd_tab, king, file = rdat)
} else {
  save(dat.inter, rel_tab, fam_table, ibd_tab, king, file = rdat)
}
