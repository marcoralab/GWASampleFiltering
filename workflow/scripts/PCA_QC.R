#!/usr/bin/env Rscript

message("Loading packages")

suppressPackageStartupMessages(library(dplyr))
library(readr)
library(magrittr)
suppressPackageStartupMessages(library(tidyr))
library(stringr)
library(tibble)
suppressPackageStartupMessages(library(purrr))

# Get geometric median
## rdocumentation.org/packages/bigutilsr/versions/0.3.3/topics/geometric_median

geometric_median <- function(u, tol = 1e-10, maxiter = 1000, by_grp = NULL) {
  if (!is.null(by_grp))
    return(do.call("rbind", by(u, by_grp, geometric_median)))
  u_old <- colMeans(u)
  for (k in seq_len(maxiter)) {
    norm <- sqrt(rowSums(sweep(u, 2, u_old, "-")^2))
    u_new <- colSums(sweep(u, 1, norm, "/")) / sum(1 / norm)
    diff <- max(abs(u_new - u_old))
    if (diff < tol)
      break
    u_old <- u_new
  }
  if (k == maxiter)
    warning("The maximum number of iterations has been reached.")
  u_new
}

# assign sample to cluster
## https://www.biorxiv.org/content/10.1101/2020.10.06.328203v2.full
## adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list

find_cluster <- function(df, clusters) {
  superpops <- clusters$superpop
  samp_pcs <- select(df, starts_with("PC"))
  mat <- bind_rows(clusters, samp_pcs) %>% {suppressWarnings(dist(.))}
  # mat
  clus <- which.min(as.matrix(mat)[6, 1:5])
  dplyr::mutate(df, superpop_infered = superpops[clus])
}

vec <- snakemake@input[["eigenvec"]]
val <- snakemake@input[["eigenval"]]
base <- snakemake@input[["pops"]]
target <- snakemake@input[["fam"]]
output <- snakemake@output[["excl"]]
rmd <- snakemake@output[["rmd"]]
sample <- snakemake@wildcards[["sample"]]
population <- snakemake@params[["superpop"]]
extraref <- snakemake@params[["extraref"]]
sdev <- snakemake@params[["sd"]]
pcs_out_path <- snakemake@output[["pcs_pops"]]
save.image(file = paste0(output, ".params.Rdata"))

#load("output/ADGC/x_present_AA/ADC8-AA_exclude.pca.params.Rdata")

if (tolower(population) == "all") {
  all_pops <- T
} else {
  all_pops <- F
}

##---- Read in Data ----##
message("Reading data files")

# count columns and PCs
n_eig <- count_fields(vec, tokenizer_delim(delim = " "), n_max = 1) - 2

# Generate colnames
pc_names <- paste0("PC", 1:n_eig)
names_col <- c("FID", "IID", pc_names)

# Read in eigenvectors and z-score transform
pca_orig <- read_delim(vec,
                  delim = " ", col_names = names_col,
                  col_types = cols(.default = "d", FID = "c", IID = "c")) %>%
         mutate_at(pc_names, function(x) as.vector(scale(x)))

# read in egienvalues
eigenval <- val %>%
  read_lines %>%
  parse_number %>%
  tibble(eigenval = .,
         PC = factor(pc_names, levels = pc_names)) %>% #PC Names
  mutate(PVE = round(eigenval / sum(eigenval), 3)) %>% #PVE
  select(PC, eigenval, PVE) #Reorder columns

# population data file, usually from 1000 genomes and potentially with extra ref
base_pops_raw <- read_table2(base, col_types = cols(.default = "c"))

# population data from target set
famcols <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")
target_pops_raw <- read_table2(target, col_names = famcols,
  col_types = "ccccii")

message("Processing data")
# ---- Data wrangling ---- #

# Read in populations and superpops
tg_pops <- read_tsv("workflow/resources/tg_subpops.tsv", col_types = "cccc")
populations <- tg_pops %>% select(pop, spop) %>% deframe %>% as.list
superpops <- unlist(populations) %>% unique()

extra_pops <- base_pops_raw[
    !(base_pops_raw$Population %in% names(populations)), ] %>%
  distinct(Population) %>%
  pull(Population)

if (length(extra_pops) == 1 && extraref == extra_pops) {
  populations[extraref] <- population
} else if (length(extra_pops) == 0 && extraref != "none") {
  stop("Extra population missing from reference!")
} else if (length(extra_pops) > 1) {
  stop("Non-1kG population codes are not yet implemented. Go bug Brian.")
}

# Deal with invalid cohort names

if (sample %in% names(populations)) {
  sample <- paste0("s_", sample)
}

if (sample %in% populations) {
  sample_s <- paste0("s_", sample)
} else {
  sample_s <- sample
}


##  Munge population dataframes from 1000 genomes
base_pops <- base_pops_raw %>%
  mutate(cohort = "Reference",
         superpop = recode(.$Population, !!!populations))

##  Munge target population dataframes
target_pops <- target_pops_raw %>%
  select(FID, IID) %>%
  mutate(Population = sample, superpop = sample_s,
    cohort = sample_s)

## Check this
remove_tg <- TRUE
if (remove_tg) {
  target_pops <- target_pops %>%
    filter(!(IID %in% base_pops$IID & FID %in% base_pops$FID))
}

# fix improperly split FID_IID
pca_fidiid <- pca_orig %>%
  unite("FIDIID", FID, IID, sep = "_")


##  Munge PCA, base pop and target pop
both_pops <- target_pops %>%
  bind_rows(base_pops) %>%
  ##### FIX BAD FID_IID SPLIT #####
  unite("FIDIID", FID, IID, sep = "_", remove = F)

pca_corrected <- pca_fidiid %>%
  left_join(both_pops, by = "FIDIID") %>%
  select(any_of(names(both_pops)), everything(), -FIDIID) %>%
  mutate(FID = str_remove(FID, "^1000g___"))

## Colours for plots
pca_col <- pca_corrected %>%
  count(superpop) %>%
  mutate(color = ifelse(superpop == sample_s, "black", NA)) %>%
  mutate(color = ifelse(superpop == "AFR", "#E69F00", color)) %>%
  mutate(color = ifelse(superpop == "AMR", "#0072B2", color)) %>%
  mutate(color = ifelse(superpop == "EAS", "#009E73", color)) %>%
  mutate(color = ifelse(superpop == "EUR", "#CC79A7", color)) %>%
  mutate(color = ifelse(superpop == "SAS", "#D55E00", color)) %>%
  mutate(color = ifelse(superpop == "MID", "#56B4E9", color)) %>%
  mutate(color = ifelse(superpop == "SAS", "#F0E442", color))

# ternary plot and assignment:

message("Ternary Plots")

# Pull out 1000 genomes samples
kg <- filter(pca_corrected, cohort == "Reference")

# find geometric median of each PC for each cluster
clusters <-
  select(kg, starts_with("PC")) %>%
  geometric_median(by_grp = kg$superpop) %>%
  as_tibble(rownames = "superpop")

# extract sample information and assign to cluster
pca <- pca_corrected %>%
  group_split(IID) %>%
  map_df(find_cluster, clusters)

write_tsv(pca, pcs_out_path)

report_settings <- list(
  all_pops = all_pops,
  filter_inference = F,
  filter_sd = F
)

if (all_pops) {
  tab_1 <- as.data.frame(clusters)
  pca_sample <- pca %>%
    filter(cohort != "Reference") %>%
    mutate(pop_outliers = F)
  tab_pop_exclusions <- tibble(Exclusions = "None")
} else if (!is.numeric(sdev)) {
  report_settings$filter_inference <- T
  tab_1 <- as.data.frame(clusters)
  pca_sample <- pca %>%
    filter(cohort != "Reference") %>%
    mutate(pop_outliers = !(superpop_infered %in% population))
  tab_pop_exclusions <- pca_sample %>%
    filter(pop_outliers) %>%
    select(FID, IID, "Infered Superpopulation" = superpop_infered)
} else {
  report_settings$filter_sd <- T
  # ---- Population Outliers ---- #
  # Calculate the mean and ± specified number of sd for each PC of ref pop

  message("Table 1")

  chosen_pca <- pca_corrected %>%
    gather(key = "PC", value = "eigenvalue", !!paste0("PC", 1:n_eig)) %>%
    filter(superpop == population) %>%
    group_by(superpop, PC) %>%
    dplyr::summarize(mean = mean(eigenvalue), sd = sd(eigenvalue),
                     .groups = "drop") %>%
    mutate(lower = mean - sd * sdev, upper = mean + sd * sdev) %>%
    mutate(PC = factor(PC, levels = paste0("PC", 1:n_eig))) %>%
    arrange(PC)

  #***Table 1:*** Mean and SD of PC in chosen population
  tab_1 <- as.data.frame(chosen_pca)

  # For each sample individual, determine if ± specified SD from chosen pop
  # for each principal component

  message("Table 2")

  pca_range <- function(pc, vals) {
    lower <- chosen_pca[chosen_pca$PC == pc, "lower"] %>% unlist %>% unname
    upper <- chosen_pca[chosen_pca$PC == pc, "upper"] %>% unlist %>% unname
    sapply(vals, function(x) x < lower | x > upper) %>% as.logical
  }

  mut_pca <- function(df, pc, pc_use) {
    pc <- 1:pc
    pc_use <- paste0("PC", pc, ".outliers")[1:pc_use]
    for (i in paste0("PC", pc)) {
      df %<>% mutate(`!!`(paste0(i, ".outliers")) := pca_range(i, .[, i]))
    }
    dplyr::mutate(df, pop_outliers = rowSums(dplyr::select(df, !!pc_use)) > 0)
  }

  pca_sample <- pca_corrected %>%
    filter(superpop == sample_s) %>%
    mut_pca(n_eig, 10)

  # ***Table 2:*** Population Outliers for each PC

  outlier_cols <- paste0("PC", 1:n_eig, ".outliers")
  outlier_cols_rename <- outlier_cols
  names(outlier_cols_rename) <- paste0("PC", 1:n_eig)

  replace_tf <- function(x) ifelse(x == T, "Yes", ifelse(x == F, "No", x))

  tab_pop_exclusions <- pca_sample %>%
    select(FID, IID, !!outlier_cols, pop_outliers) %>%
    filter(pop_outliers == T) %>%
    dplyr::rename(!!!outlier_cols_rename) %>%
    select(-pop_outliers) %>%
    transmute_all(replace_tf)
}

# ---- Write out outliers ---- #

no_outliers <- sum(pca_sample$pop_outliers) == 0
exclude_pop_outliers <- pca_sample %>%
  filter(pop_outliers == TRUE) %>%
  select(FID, IID)

write_tsv(exclude_pop_outliers, output, col_names = F)

save(tab_1, tab_pop_exclusions, no_outliers, pca, pca_sample, eigenval,
     pca_col, sdev, report_settings, file = rmd)
