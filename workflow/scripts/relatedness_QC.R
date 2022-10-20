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
outfile_allrelate <- snakemake@output[["allrelate"]]
rdat <- snakemake@output[["rdat"]]
# save.image(file = paste0(rdat, '.rda'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Read in Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

norel <- any(read_lines(snakemake@input[["genome"]]) == "norel")

fam <- snakemake@input[["fam"]] %>%
  read_table(col_types = "cc---i", col_names = c("FID", "IID", "status"))

if (snakemake@input[["exclude"]] != "/dev/null") {
  exclusions <- snakemake@input[["exclude"]] %>%
    read_tsv(col_types = "ccc") %>%
    group_by(FID, IID) %>%
    # collapse all exclusions into one row and separate the reasons with commas
    dplyr::summarise(why = paste(why, collapse = ","), .groups = "drop") %>%
    mutate(qc_failed = TRUE)
  fam <- fam %>%
    left_join(exclusions, by = c("FID", "IID")) %>%
    mutate(qc_failed = replace_na(qc_failed, FALSE))
} else {
  fam <- mutate(fam, qc_failed = FALSE)
}

dat_inter_all_kin0 <- paste0(genome_file, ".all.kin0")
dat_inter_all_kin <- paste0(genome_file, ".all.kin")
dat_inter_kin0_path <- paste0(genome_file, ".kin0")
dat_inter_kin_path <- paste0(genome_file, ".kin")

kin0 <- file.exists(dat_inter_all_kin0)
kin <- file.exists(dat_inter_all_kin)
kin0_rel_exists <- file.exists(dat_inter_kin0_path)
kin_rel_exists <- file.exists(dat_inter_kin_path)

if (kin0) {
  exclude_unrelated <- file.size(dat_inter_all_kin0) >= 4e+10
} else {
  exclude_unrelated <- FALSE
}

if (exclude_unrelated) {
  read_kinship <- function(file, ...) {
    filter_unrelate <- function(x, i) {
      filter(x, Kinship > 0)
    }
    read_tsv_chunked(file, DataFrameCallback$new(filter_unrelate),
                     chunk_size = 100000, ...)
  }
} else {
  read_kinship <- read_table
}

if (kin0) {
  dat_inter_all_kin0 <- dat_inter_all_kin0 %>%
    read_kinship(col_types = cols(
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

if (!norel && kin0_rel_exists) {
  dat_inter_kin0 <- dat_inter_kin0_path %>%
    read_kinship(col_types = cols(
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

if (kin) {
  dat_inter_all_kin <- dat_inter_all_kin %>%
    read_kinship(col_types = cols(
      .default = col_double(),
      FID = col_character(),
      ID1 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
    )) %>%
    rename(IID1 = ID1, IID2 = ID2, FID1 = FID) %>%
    mutate(FID2 = FID1, PI_HAT = ifelse(Kinship > 0, 2 * Kinship, 0))
}

if (!norel && kin_rel_exists) {
  dat_inter_kin <- dat_inter_kin_path %>%
    read_kinship(col_types = cols(
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

if (!norel && kin0_rel_exists && kin_rel_exists) {
  print("Relatedness kin and kin0 files present")
  dat_inter <- bind_rows(dat_inter_kin0, dat_inter_kin) %>%
    distinct(FID1, IID1, FID2, IID2, .keep_all = TRUE)
} else if (!norel && kin0_rel_exists) {
  print("Relatedness kin0 file present")
  dat_inter <- dat_inter_kin0
} else if (!norel && kin_rel_exists) {
  print("Relatedness kin file present")
  dat_inter <- dat_inter_kin
}

if (kin0 && kin) {
  print("Kin and kin0 files present")
  dat_inter_all <- bind_rows(dat_inter_all_kin0, dat_inter_all_kin) %>%
    distinct(FID1, IID1, FID2, IID2, .keep_all = TRUE)
} else if (kin0) {
  print("Kin0 file present")
  dat_inter_all <- dat_inter_all_kin0
} else if (kin) {
  print("Kin file present")
  dat_inter_all <- dat_inter_all_kin
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
    mutate(
      relationship = case_when(
        InfType %in% c("1st", "2nd", "3rd", "4th") ~
          paste(InfType, "degree"),
        InfType == "UN" ~ "unrelated",
        InfType == "Dup/MZ" ~ "identical-twins",
        InfType == "PO" ~ "parent-child",
        InfType == "FS" ~ "full-siblings",
        TRUE ~ InfType))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  Exclude Samples
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# select samples with kinship cofficents > 0.1875
# https://link.springer.com/protocol/10.1007/978-1-60327-367-1_19

remove_samples <- function(ibdcoeff, fam, msg = "closely related to") {
  fam_fi <- fam %>%
    mutate(FI = paste0(FID, "_-_-tempsep-_-_", IID)) %>%
    mutate(status = ifelse(status > 2, 0.5, status))

  ibdcoeff %<>%
    mutate(FI1 = paste0(FID1, "_-_-tempsep-_-_", IID1),
           FI2 = paste0(FID2, "_-_-tempsep-_-_", IID2))
  related_samples <- NULL
  excluded <- c()
  fam_table <- tibble(FID = c("deleteme"),
                      IID = c("deleteme"),
                      Related = c("deleteme"))
  while (nrow(ibdcoeff) > 0) {
    test_tab <- plyr:::count(c(ibdcoeff$FI1, ibdcoeff$FI2))
    if (!("x" %in% names(test_tab))) {
      print(ibdcoeff)
    }
    sample.counts <- plyr:::count(c(ibdcoeff$FI1, ibdcoeff$FI2)) %>%
      as_tibble %>%
      rename(FI = x) %>%
      mutate(FI = as.character(FI)) %>%
      inner_join(fam_fi, by = "FI") %>%
      arrange(desc(qc_failed), status, desc(freq))
    rm.sample <- sample.counts[[1, "FI"]]
    id_ <- str_split(rm.sample, "_-_-tempsep-_-_")[[1]]
    fid <- id_[1]
    iid <- id_[2]
    remtxt <- sprintf("%s %i other samples.",
                      msg,
                      sample.counts[[1, "freq"]])
    message(paste("Removing sample", iid, remtxt))
    ft <- tibble(FID = fid, IID = iid, Related = remtxt)
    fam_table <- fam_table %>%
      bind_rows(ft)
    ibdcoeff <- ibdcoeff[ibdcoeff$FI1 != rm.sample &
                           ibdcoeff$FI2 != rm.sample, ]
    related_samples <- c(as.character(rm.sample), related_samples)
  }
  return(
    list(related_samples = related_samples,
         fam_table = filter(fam_table, Related != "deleteme"),
         exclude_samples = tibble(FI = as.character(related_samples)) %>%
           separate(FI, c("FID", "IID"), sep = "_-_-tempsep-_-_")))
}

if (norel) {
  message("No related samples.")
  fam_table <- tibble(message = "No related samples.") %>% as.data.frame()
  related_samples <- tibble(FID = character(), IID = character())
  exclude_samples <- tibble(FID = character(), IID = character())
  ibd_tab <- tibble(message = "No related samples.")
} else {
  ibdcoeff <- dat_inter %>%
    filter(PI_HAT > threshold)

  ibd_tab_create <- . %>%
    select(FID1, IID1, FID2, IID2, IBS0, Kinship, PI_HAT, relationship)

  # Iterativly remove subjects with the highest number of
  #   pairwise kinship cofficents > threshold
  # see http://www.stat-gen.org/tut/tut_preproc.html

  if (family == FALSE || family == "F") {
    ibd_tab <- ibd_tab_create(ibdcoeff) %>%
      select(FID1, IID1, FID2, IID2, IBS0, Kinship, PI_HAT, relationship)
    message("Working with unrelated samples.")
    removal <- remove_samples(ibdcoeff, fam)
  } else {
    ibdcoeff_dup <- dat_inter %>%
      filter(relationship == "identical-twins")
    ibd_tab <- ibd_tab_create(ibdcoeff_dup)
    message("Working with related samples.")
    message("Finding duplicates/twins")
    removal <- remove_samples(ibdcoeff_dup, fam, "duplicate/twin of")
    message("\nFinding all relationships")
    removal_allrelate <- remove_samples(ibdcoeff, fam)
  }
  # Don't simplify fam_table. It will break.
  fam_table <- removal$fam_table
  exclude_samples <- removal$exclude_samples
}

dat_forplots <- dat_inter_all

##  write out samples to be excluded
write_tsv(exclude_samples, outfile, col_names = TRUE)
if (norel) {
  save(norel, dat_forplots, rel_tab, fam_table, ibd_tab, file = rdat)
} else if (family == TRUE || family == "T") {
  write_tsv(removal_allrelate$exclude_samples,
            outfile_allrelate, col_names = TRUE)
  save(norel, dat_inter, dat_forplots, rel_tab, fam_table, ibd_tab,
       exclude_unrelated, removal_allrelate, file = rdat)
} else {
  save(norel, dat_inter, dat_forplots, rel_tab, fam_table, ibd_tab,
       exclude_unrelated, file = rdat)
}
