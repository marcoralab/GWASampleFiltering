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
  "stringr"
  ))

##--- Set Params ---##
# parser <- ArgumentParser()
# # specify our desired options
# # by default ArgumentParser will add an help option
# parser$add_argument("--vec", help = "Eigenvectors file", required = T)
# parser$add_argument("--val", help = "Eigenvalues file", required = T)
# parser$add_argument("-b", "--base", help = "Ped file from 1kgp", required = T)
# parser$add_argument("-t", "--target",
#   help = "Fam file for sample. Can include 1kgp.", required = T)
# parser$add_argument("--TGremove", action="store_true")
# parser$add_argument("-s", "--sample", help = "Sample name", required = T)
# parser$add_argument("-o", "--output", help = "Sample exclusion file")
# parser$add_argument("-p", "--population", default = "EUR",
#   help = "Superpop [Default: EUR; Options: EUR, AMR, AFR, EAS, SAS]")
# parser$add_argument("-R", "--rmd", help = "Rdata file for RMD report")
#
# params <- parser$parse_args()
# save(params, file=paste0(params$output, ".params.Rdata"))

vec = snakemake@input[['eigenvec']]
val = snakemake@input[['eigenval']]
base = snakemake@input[['tgped']]
target = snakemake@input[['fam']]
output = snakemake@output[['excl']]
rmd = snakemake@output[['rmd']]
sample = snakemake@params[['samp']]
population = snakemake@params[['superpop']]

##  Standarize varibles to have a mean of 0 and sd of 1
zscore = function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)}

##---- Read in Data ----##
message("Reading data files")

# count colums and PCs
n_cols <- count_fields(vec, tokenizer_delim(delim = " "), n_max = 1)
n_eig <- n_cols - 2

# Generate colnames
pc_names <- paste0("PC", 1:n_eig)
names_col <- c("FID", "IID", pc_names)

# Read in eigenvectors and z-score transform
pca_orig <- read_delim(vec,
                  delim = " ", col_names = names_col,
                  col_types = cols(.default = "d", FID = "c", IID = "c")) %>%
         mutate_at(pc_names, zscore)

# read in egienvalues
eigenval.raw <- parse_number(read_lines(val))
eigenval <- tibble(eigenval = eigenval.raw,
                   PC = factor(pc_names, levels = pc_names)) %>% #PC Names
              mutate(PVE = round(eigenval / sum(eigenval), 3)) %>% #PVE
              select(PC, eigenval, PVE) #Reorder columns

# population data file from 1000 genomes
ct <- cols_only(`Individual ID` = "c",
                `Family ID` = "c",
                Population = "c")
base_pops.raw <- read_tsv(base, col_types = ct)

# population data from target set
famcols = c("FID", "IID", "PID", "MID", "Sex", "Pheno")
target_pops.raw <- read_table2(target, col_names = famcols,
  col_types = "ccccii")

message("Processing data")
# ---- Data wrangling ---- #

s_pops <- list(GBR = "EUR", FIN = "EUR", IBS = "EUR", CEU = "EUR", TSI = "EUR",
               CHS = "EAS", CDX = "EAS", KHV = "EAS", CHB = "EAS", JPT = "EAS",
               PJL = "SAS", BEB = "SAS", GIH = "SAS", STU = "SAS", ITU = "SAS",
               PUR = "AMR", PEL = "AMR", CLM = "AMR", MXL = "AMR",
               ACB = "AFR", GWD = "AFR", ESN = "AFR", MSL = "AFR", LWK = "AFR",
               YRI = "AFR", ASW = "AFR")



# Deal with invalid cohort names

if (sample %in% names(s_pops)) {
  sample <- paste0("s_", sample)
}

if (sample %in% s_pops) {
  sample_s <- paste0("s_", sample)
} else {
  sample_s <- sample
}


##  Munge population dataframes from 1000 genomes
base_pops <- base_pops.raw %>%
  dplyr::rename(FID = `Family ID`, IID = `Individual ID`) %>%
  mutate(cohort = "1kgenomes",
         superpop = recode(.$Population, !!!s_pops))

##  Munge target population dataframes
target_pops <- target_pops.raw %>%
  select(FID, IID) %>%
  mutate(Population = sample, superpop = sample_s,
    cohort = sample_s)

## Check this
TGremove = TRUE
if ( TGremove ) {
  target_pops <- target_pops %>%
    filter(!(IID %in% base_pops$IID))
}

# fix improperly split FID_IID
pca_FIDIID <- pca_orig %>%
  unite("FIDIID", FID, IID, sep = "_")


##  Munge PCA, base pop and target pop
pca <- target_pops %>%
  bind_rows(base_pops) %>%
  ##### FIX BAD FID_IID SPLIT #####
  unite("FIDIID", FID, IID, sep = "_", remove = F) %>%
  right_join(pca_FIDIID, by = "FIDIID") %>%
  select(-FIDIID) %>%
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
# Calculate the mean and ± 6 sd for each PC of 1000 Genomes EUR population

chosen_pca <- pca %>%
  gather(key = "PC", value = "eigenvalue", !!paste0("PC", 1:n_eig)) %>%
  filter(superpop == population) %>%
  group_by(superpop, PC) %>%
  dplyr::summarize(mean = mean(eigenvalue), sd = sd(eigenvalue)) %>%
  mutate(lower = mean - sd * 6, upper = mean + sd*6) %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:n_eig))) %>%
  arrange(PC)

#***Table 1:*** Mean and SD of PC in chosen population
tab_1 <- as.data.frame(chosen_pca)

# For each individual in the sample, determine if ±6SD from the chosen population
# for each principal component

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

#	---- Write out outliers ---- #
exclude_pop_outliers <- sample.pca %>%
  filter(pop.outliers == TRUE) %>%
  select(FID, IID)

write_tsv(exclude_pop_outliers, output, col_names = F)
save(tab_1, tab_2, no_outliers, pca, sample.pca, eigenval, pca_col,
     file = rmd)
