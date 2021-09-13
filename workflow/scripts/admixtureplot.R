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
  "\noutput:" ,p.out,"\n"
)
message ("Loading packages")
library(dplyr)
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

Q.dat.long <- Q.dat %>%
  rowwise(IID) %>%
  mutate(
    maxval=max(c_across(c(AFR, SAS, EAS, EUR, AMR))),
    matchval=which.max(c_across(c(AFR, SAS, EAS, EUR, AMR))),
  ) %>%
  ungroup() %>%
  mutate(corder = order(matchval, -maxval)) %>%
  pivot_longer(c("AFR", "SAS", "EAS", "EUR", "AMR"), names_to = "K", values_to = "prop" ) %>%
  arrange(matchval, -maxval) %>%
  mutate(IID = fct_inorder(IID))

colour <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
plotQ <- ggplot(Q.dat.long , aes(x = IID, y = prop, fill = K)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  scale_fill_manual(name="Super Population",
                    values=colour ,
                    breaks=c("AFR","AMR","EAS","EUR","SAS") ,
                    labels=c("AFR", "AMR","EAS", "EUR", "SAS")) +
  theme_classic() + labs(x = "Indivuals", y = "Global Ancestry", color ="Super Population") +
  facet_grid(~fct_relevel(pca_super_pop, "EUR", "AFR", "AMR", "EAS", "SAS"), switch = "x", scales = "free", space = "free")+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x =element_blank(),
    panel.grid.major.x = element_blank())


ggsave(p.out, plot = plotQ, width = 12, height = 6, units = "in")
