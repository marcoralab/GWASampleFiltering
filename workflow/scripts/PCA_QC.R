#!/usr/bin/env Rscript
## ---- Load Required R packages ---- ##
import <- Vectorize(function(x) {
  suppressWarnings(suppressMessages(require(x, character.only = T, quietly = T,
    warn.conflicts = F)))
})

message("Loading packages")

status <- import(c(
  "dplyr",
  "readr",
  "magrittr",
  "argparse",
  "tidyr",
  "stringr",
  "tibble",
  "purrr"
  ))

`%nin%` <- Negate(`%in%`)

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
## http://adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list/
find_cluster <- function(df, clusters) {
  superpops <- clusters$superpop
  samp_pcs <- select(df, starts_with("PC"))
  mat <- bind_rows(clusters, samp_pcs) %>% {suppressWarnings(dist(.))}
  # mat
  clus <- as.matrix(mat)[6, 1:5] %>% which.min()
  df %>% mutate(superpop_infered = superpops[clus])
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

#load("/sc/arion/projects/LOAD/Brian/src/GWASampleFiltering/output/ADGC/x_present_AA/ADC8-AA_exclude.pca.params.Rdata")

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
pca_FIDIID <- pca_orig %>%
  unite("FIDIID", FID, IID, sep = "_")


##  Munge PCA, base pop and target pop
both_pops <- target_pops %>%
  bind_rows(base_pops) %>%
  ##### FIX BAD FID_IID SPLIT #####
  unite("FIDIID", FID, IID, sep = "_", remove = F)

pca <- pca_FIDIID %>%
  left_join(both_pops, by = "FIDIID") %>%
  select(any_of(names(both_pops)), everything(), -FIDIID) %>%
  mutate(FID = str_remove(FID, "^1000g___"))

## Colours for plots
pca_col <- pca %>%
  count(superpop) %>%
  mutate(colour = ifelse(superpop == sample_s, "black", NA)) %>%
  mutate(colour = ifelse(superpop == "AFR", "#E41A1C", colour)) %>%
  mutate(colour = ifelse(superpop == "AMR", "#377EB8", colour)) %>%
  mutate(colour = ifelse(superpop == "EAS", "#4DAF4A", colour)) %>%
  mutate(colour = ifelse(superpop == "EUR", "#984EA3", colour)) %>%
  mutate(colour = ifelse(superpop == "SAS", "#FF7F00", colour))

# ---- Population Outliers ---- #
# Calculate the mean and ± specified number of sd for each PC of 1000 Genomes EUR population

message("Table 1")

chosen_pca <- pca %>%
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

# For each individual in the sample, determine if ± specified SD from the chosen population
# for each principal component

message("Table 2")

pca_range <- function(PC, vals) {
  lower <- chosen_pca[chosen_pca$PC == PC, "lower"] %>% unlist %>% unname
  upper <- chosen_pca[chosen_pca$PC == PC, "upper"] %>% unlist %>% unname
  sapply( vals, function(x) {x < lower | x > upper} ) %>% as.logical
}

mut_pca <- function(df, PC, PC_use) {
  PC <- 1:PC
  Po_use <- paste0("PC", PC, ".outliers")[1:PC_use]
  for (i in paste0("PC", PC)) {
    df %<>% mutate(`!!`(paste0(i, ".outliers")) := pca_range(i, .[, i]))
  }
  df %>% dplyr::mutate(pop.outliers = rowSums(dplyr::select(., !!Po_use)) > 0)
}

sample.pca <- pca %>%
  filter(superpop == sample_s) %>%
  mut_pca(n_eig, 10)

n_outliers <- sum(sample.pca$pop.outliers)
no_outliers <- n_outliers == 0

# ***Table 2:*** Population Outliers for each PC

outlier_cols <- paste0("PC", 1:n_eig, ".outliers")
outlier_cols_rename <- outlier_cols
names(outlier_cols_rename) <- paste0("PC", 1:n_eig)

tab_2 <- as.data.frame(
  sample.pca %>%
  select(FID, IID, !!outlier_cols, pop.outliers) %>%
  filter(pop.outliers == T) %>%
  dplyr::rename(!!!outlier_cols_rename, Outlier = pop.outliers)
)

# ternary plot and assignment:

message("Ternary Plots")

pca_geomed <- pca_orig %>%
  left_join(both_pops, by = c("FID", "IID")) %>%
  mutate(superpop = ifelse(is.na(superpop), "sample", superpop),
       Population = ifelse(is.na(Population), "sample", Population))

# Pull out 1000 genomes samples
kg <- filter(pca_geomed, cohort == "Reference")

# find geometric median of each PC for each cluster
clusters <-
  select(kg, starts_with("PC")) %>%
  geometric_median(by_grp = kg$superpop) %>%
  as_tibble(rownames = "superpop")

# extract sample information and assign to cluster
pcs_geomed <- pca_geomed %>%
  group_split(IID) %>%
  map_df(find_cluster, clusters)

write_tsv(pcs_geomed, pcs_out_path)

# ---- Write out outliers ---- #
exclude_pop_outliers <- sample.pca %>%
  filter(pop.outliers == TRUE) %>%
  select(FID, IID)

write_tsv(exclude_pop_outliers, output, col_names = F)
save(tab_1, tab_2, no_outliers, pca, sample.pca, eigenval, pca_col, sdev,
     pcs_geomed, file = rmd)
