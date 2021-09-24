#!/usr/bin/env Rscript
## ========================================================================== ##
## Make plots for model-based clustering of ancestral populations w/ ggplot
## ========================================================================== ##

## Infiles
q.path = snakemake@input[["Qdat"]]
## Outfile
p.out = snakemake@output[[1]] # Output

message(
  "\ninput: ", q.path,
  "\noutput:", p.out, "\n"
)

message ("Loading packages")
suppressPackageStartupMessages(library(dplyr))
library(tibble)
library(purrr)
library(forcats)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(readr)
library(tidyr)
library(RColorBrewer)


message("Reading admixture output file  \n")
Q.dat <- read_tsv(q.path)

spop_cols <- intersect(c("EUR", "AFR", "AMR", "EAS", "SAS"), names(Q.dat))

Q.dat.long <- Q.dat %>%
  rowwise(IID) %>%
  mutate(
    maxval=max(c_across(spop_cols)),
    matchval=which.max(c_across(spop_cols)),
  ) %>%
  ungroup() %>%
  mutate(corder = order(matchval, -maxval)) %>%
  pivot_longer(spop_cols, names_to = "K", values_to = "prop" ) %>%
  arrange(matchval, -maxval) %>%
  mutate(IID = fct_inorder(IID))

colour <- c(AFR = "#E41A1C", AMR = "#377EB8", EAS = "#4DAF4A",
            EUR = "#984EA3", SAS = "#FF7F00")[spop_cols]
plotQ <- ggplot(Q.dat.long , aes(x = IID, y = prop, fill = K)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  scale_fill_manual(name="Super Population",
                    values=colour) +
  theme_classic() + labs(x = "Indivuals", y = "Global Ancestry", color ="Super Population") +
  facet_grid(~fct_relevel(pca_super_pop, spop_cols), switch = "x", scales = "free", space = "free")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x =element_blank(),
    panel.grid.major.x = element_blank())


ggsave(p.out, plot = plotQ, width = 12, height = 6, units = "in")
