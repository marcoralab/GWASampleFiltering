library(tidyverse)
`%nin%` = negate(`%in%`)

# Get geometric median
## https://www.rdocumentation.org/packages/bigutilsr/versions/0.3.3/topics/geometric_median
geometric_median <- function (U, tol = 1e-10, maxiter = 1000, by_grp = NULL){
  if (!is.null(by_grp))
    return(do.call("rbind", by(U, by_grp, geometric_median)))
  u.old <- colMeans(U)
  for (k in seq_len(maxiter)) {
    norm <- sqrt(rowSums(sweep(U, 2, u.old, "-")^2))
    u.new <- colSums(sweep(U, 1, norm, "/"))/sum(1/norm)
    diff <- max(abs(u.new - u.old))
    if (diff < tol)
      break
    u.old <- u.new
  }
  if (k == maxiter)
    warning("The maximum number of iterations has been reached.")
  u.new
}

# assign sample to cluster
## https://www.biorxiv.org/content/10.1101/2020.10.06.328203v2.full
## http://adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list/
find_cluster <- function(df){
  iid <- select(df, starts_with("PC"))
  mat <- bind_rows(clusters, iid) %>% dist(.)
  # mat
  clus <- as.matrix(mat)[6,1:5] %>% which.min()
  df %>% mutate(super_pop = pops[clus])
}

# ref_pops.path <- "GWASampleFiltering/reference/1kG_pops.txt"
# ref_superpops.path <- "GWASampleFiltering/reference/1kG_superpops.txt"
# vec.path <- "GWASampleFiltering/sandbox/output/gsa1234567_1kG_merged.eigenvec"

ref_pops.path <- snakemake@input[['ref_pops']]
vec.path <- snakemake@input[['eigenvec']]
pcs.out.path <- snakemake@output[['pcs_pops']]

# population names
s_pops <- list(GBR = "EUR", FIN = "EUR", IBS = "EUR", CEU = "EUR", TSI = "EUR",
               CHS = "EAS", CDX = "EAS", KHV = "EAS", CHB = "EAS", JPT = "EAS",
               PJL = "SAS", BEB = "SAS", GIH = "SAS", STU = "SAS", ITU = "SAS",
               PUR = "AMR", PEL = "AMR", CLM = "AMR", MXL = "AMR",
               ACB = "AFR", GWD = "AFR", ESN = "AFR", MSL = "AFR", LWK = "AFR",
               YRI = "AFR", ASW = "AFR")

pops <- unlist(s_pops) %>% unique()

ref <- read_table2(ref_pops.path) %>%
  mutate(cohort = "Reference",
         superpop = recode(.$Population, !!!s_pops))

# Formating data
vec <- read_table2(vec.path, col_names = F) %>%
  rename(fid = X1, iid = X2, PC1 = X3, PC2 = X4, PC3 = X5, PC4 = X6, PC5 = X7, PC6 = X8, PC7 = X9, PC8 = X10, PC9 = X11, PC10 = X12) %>%
  left_join(ref, by = c('iid' = 'IID')) %>%
  mutate(superpop = ifelse(is.na(superpop), "sample", superpop),
         Population = ifelse(is.na(Population), "sample", Population))

# Pull out 1000 genomes samples
kg <- filter(vec, iid %in% ref$IID)

# find geometric median of each PC for each cluster
clusters <- kg %>%
  group_split(superpop) %>%
  magrittr::set_names(pops) %>%
  map(., select, starts_with("PC")) %>%
  map(., geometric_median) %>%
  bind_rows(.id = "iid")

# extract sample information and assign to cluster
samples <- vec %>%
  filter(., iid %nin% ref$IID) %>%
  group_split(iid) %>%
  map_df(., find_cluster)

pcs <- bind_rows(kg, samples)

write_tsv(pcs, pcs.out.path)
