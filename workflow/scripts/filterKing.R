#!/usr/bin/env Rscript

require(tibble)
require(readr)
suppressMessages(require(dplyr))

kingstem <- snakemake@params[["indat"]]
exclude <- snakemake@input[["exclude"]]
outfile <- snakemake@output[[1]]

message(
  "kingstem: ", kingstem, "\n",
  "exclude: ", exclude, "\n",
  "outfile: ", outfile, "\n"
)

excluded <- read_tsv(exclude, col_names = c("FID", "IID"), col_types = "cc")

kin0file <- paste0(kingstem, ".kin0")
kinfile <- paste0(kingstem, ".kin")

if (nrow(excluded) == 0) {
  message("no kinship exclusions.")
  if (file.exists(kinfile)) {
    message("Symlinking kinfile")
    kinout <- paste0(kingstem, ".popfilt.kin")
    system(sprintf("ln -srf %s.kin %s", kingstem, kinout))
  } else {
    message(sprintf("%s does not exist.", kin0file))
    kinout <- ""
  }

  if (file.exists(kin0file)) {
    message("Symlinking kin0file")
    kin0out <- paste0(kingstem, ".popfilt.kin0")
    system(sprintf("ln -srf %s.kin0 %s", kingstem, kin0out))
  } else {
    message(sprintf("%s does not exist.", kin0file))
    kin0out <- ""
  }
} else {
  if (file.exists(kinfile)) {
    message(sprintf("Loading %s", kinfile))
    kincols <- cols(
      .default = col_double(),
      FID = col_character(),
      ID1 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
    )

    kinout <- paste0(kingstem, ".popfilt.kin")
    read_tsv(kinfile, col_types = kincols) %>%
      anti_join(excluded, by = c("FID", "ID1" = "IID")) %>%
      anti_join(excluded, by = c("FID", "ID2" = "IID")) %>%
      write_tsv(kinout)
  } else {
    message(sprintf("%s does not exist.", kin0file))
    kinout <- ""
  }

  if (file.exists(kin0file)) {
    message(sprintf("Loading %s", kin0file))
    kin0cols <- cols(
      .default = col_double(),
      FID1 = col_character(),
      ID1 = col_character(),
      FID2 = col_character(),
      ID2 = col_character(),
      N_SNP = col_integer()
    )

    kin0out <- paste0(kingstem, ".popfilt.kin0")
    read_tsv(kin0file, col_types = kin0cols) %>%
      anti_join(excluded, by = c("FID1" = "FID", "ID1" = "IID")) %>%
      anti_join(excluded, by = c("FID2" = "FID", "ID2" = "IID")) %>%
      write_tsv(kin0out)
  } else {
    message(sprintf("%s does not exist.", kin0file))
    kin0out <- ""
  }
}

message("Exporting... ", outfile)

if (nchar(paste0(kinout, kin0out)) > 0) {
  fileconn <- file(outfile)
  writeLines(c(kinout, kin0out), fileconn)
  close(fileconn)
} else {
  stop("No KING output files!")
}
