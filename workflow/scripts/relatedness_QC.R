#!/usr/bin/env Rscript
##	R script for identifing samples that are crypticaly related from PLINK files.

suppressMessages(require(readr))
suppressMessages(require(dplyr))
suppressMessages(require(tibble))
suppressMessages(require(stringr))
suppressMessages(require(tidyr))
suppressMessages(require(magrittr))

`%nin%` <- Negate(`%in%`)

# .genome output file from PLINK
genome_file <- snakemake@params[["geno"]]
# pi_hat threshold (0.1875 recommended)
threshold <- as.numeric(snakemake@params[["threshold"]])
# Family based study, T/F
family <- as.logical(snakemake@params[["Family"]])
# .tsv file to be output
outfile <- snakemake@output[["out"]]
rdat <- snakemake@output[["rdat"]]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Read in Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

norel <- read_lines(snakemake@input[["genome"]]) == "norel"

kin0 <- F
kin <- F
dat_inter_all_kin0 <- paste0(genome_file, ".all.kin0")
dat_inter_all_kin <- paste0(genome_file, ".all.kin")
if (file.exists(dat_inter_all_kin0)) {
  kin0 <- T
  if (!norel) {
    dat_inter_kin0 <- paste0(genome_file, ".kin0") %>%
      read_table2(col_types = cols(
        .default = col_double(),
        FID1 = col_character(),
        ID1 = col_character(),
        FID2 = col_character(),
        ID2 = col_character(),
        N_SNP = col_integer(),
        InfType = col_character()
      )) %>%
      rename(IID1 = ID1, IID2 = ID2, PI_HAT = PropIBD)
  }
  dat_inter_all_kin0 <- dat_inter_all_kin0 %>%
    read_table2(col_types = cols(
      .default = col_double(),
      FID1 = col_character(),
      ID1 = col_character(),
      FID2 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
  )) %>%
    rename(IID1 = ID1, IID2 = ID2) %>%
    mutate(PI_HAT = ifelse(Kinship > 0, 2 * Kinship, 0))
}

if (file.exists(dat_inter_all_kin)) {
  kin <- T
  if (!norel) {
    dat_inter_kin <- paste0(genome_file, ".kin") %>%
      read_table2(col_types = cols(
        .default = col_double(),
        FID = col_character(),
        ID1 = col_character(),
        ID2 = col_character(),
        N_SNP = col_integer(),
        InfType = col_character()
      )) %>%
      rename(IID1 = ID1, IID2 = ID2, FID1 = FID, PI_HAT = PropIBD) %>%
      mutate(FID2 = FID1)
  }
  dat_inter_all_kin <- dat_inter_all_kin %>%
    read_table2(col_types = cols(
      .default = col_double(),
      FID = col_character(),
      ID1 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
    )) %>%
    rename(IID1 = ID1, IID2 = ID2, FID1 = FID) %>%
    mutate(FID2 = FID1, PI_HAT = ifelse(Kinship > 0, 2 * Kinship, 0))
}

if (kin0 & kin) {
  print("Kin and kin0 files present")
  dat_inter_all <- bind_rows(dat_inter_all_kin0, dat_inter_all_kin) %>%
    distinct(FID1, IID1, FID2, IID2, .keep_all = T)
  if (!norel) {
    dat_inter <- bind_rows(dat_inter_kin0, dat_inter_kin) %>%
      distinct(FID1, IID1, FID2, IID2, .keep_all = T)
  }
} else if (kin0) {
  print("Kin0 files present")
  dat_inter_all <- dat_inter_all_kin0
  if (!norel) dat_inter <- dat_inter_kin0
} else if (kin) {
  print("Kin files present")
  dat_inter_all <- dat_inter_all_kin
  if (!norel) dat_inter <- dat_inter_kin
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
                               ifelse(relationship == "first-cousin",
                                      "3rd degree", relationship)))

# Match Pi-Hat, Z0, Z1, Z2 to closet values in IBD relationship table

closest <- function(vals, ref) {
  fc <- Vectorize(function(x) {
    ref[which.min(abs(ref - x))]
  }) #finds closest
  fc(vals)
}

if (!norel) {
  dat_inter %<>%
    mutate(relationship = InfType) %>%
    mutate(relationship = ifelse(
      relationship %in% c("1st", "2nd", "3rd", "4th"),
      paste(relationship, "degree"), relationship)) %>%
    mutate(relationship = ifelse(
      relationship == "UN", "unrelated", relationship)) %>%
    mutate(relationship = ifelse(
      relationship == "Dup/MZ", "identical-twins", relationship)) %>%
    mutate(relationship = ifelse(
      relationship == "PO", "parent-child", relationship)) %>%
    mutate(relationship = ifelse(
      relationship == "FS", "full-siblings", relationship))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Exclude Samples
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# select samples with kinship cofficents > 0.1875
# https://link.springer.com/protocol/10.1007/978-1-60327-367-1_19

if (norel) {
  message("No related samples.")
  fam_table <- tibble(message = "No related samples.") %>% as.data.frame()
  related_samples <- tibble(FID = character(), IID = character())
  exclude_samples <- tibble(FID = character(), IID = character())
  ibd_tab <- tibble(message = "No related samples.")
} else {
  if (family == F | family == "F") {
    ibdcoeff <- dat_inter %>%
      filter(PI_HAT > threshold)
  } else {
    ibdcoeff <- dat_inter %>%
      filter(relationship == "identicial-twins")
  }

  ibd_tab <- ibdcoeff %>%
    select(FID1, IID1, FID2, IID2, IBS0, Kinship, PI_HAT, relationship)

  # Iterativly remove subjects with the highest number of
  #   pairwise kinship cofficents > threshold
  # see http://www.stat-gen.org/tut/tut_preproc.html

  if (family == F | family == "F") {
    message("Working with unrelated samples.")
    ibdcoeff %<>%
      mutate(FI1 = paste0(FID1, "_-_-tempsep-_-_", IID1),
             FI2 = paste0(FID2, "_-_-tempsep-_-_", IID2))
    related_samples <- NULL
    excluded <- c()
    fam_table <- tibble(FID = c("deleteme"),
                        IID = c("deleteme"),
                        Related = c("deleteme"))
    while (nrow(ibdcoeff) > 0) {
      sample.counts <- arrange(
        plyr:::count(c(ibdcoeff$FI1, ibdcoeff$FI2)), -freq)
      rm.sample <- sample.counts[1, "x"]
      id_ <- str_split(rm.sample, "_-_-tempsep-_-_")[[1]]
      fid <- id_[1]
      iid <- id_[2]
      remtxt <- sprintf("closely related to %i other samples.",
                        sample.counts[1, "freq"])
      message(paste("Removing sample", iid, remtxt))
      ft <- tibble(FID = fid, IID = iid, Related = remtxt)
      fam_table <- fam_table %>%
        bind_rows(ft)
      ibdcoeff <- ibdcoeff[ibdcoeff$FI1 != rm.sample &
                           ibdcoeff$FI2 != rm.sample, ]
      related_samples <- c(as.character(rm.sample), related_samples)
    }
    # Don't simplify fam_table. It will break.
    fam_table <- fam_table %>%
      filter(Related != "deleteme")
    exclude_samples <- tibble(FI = as.character(related_samples)) %>%
      separate(FI, c("FID", "IID"), sep = "_-_-tempsep-_-_")
  } else {
    message("Working with related samples.")
    fam_table <- as.data.frame(ibdcoeff)
    related_samples <- c(ibdcoeff$IID1, ibdcoeff$IID2)
    exclude_samples <- tibble(FID = character(), IID = character())
  }
}

##  write out samples to be excluded
write_tsv(exclude_samples, outfile, col_names = T)
if (norel) {
  save(norel, dat_inter_all, rel_tab, fam_table, ibd_tab, file = rdat)
} else {
  save(norel, dat_inter, dat_inter_all, rel_tab, fam_table, ibd_tab,
    file = rdat)
}
