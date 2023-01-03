# Split cCRE-gene pairs by chromosome

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load cCRE-gene pairs
pairs <- read_tsv(snakemake@input[[1]], col_names = FALSE, show_col_types = FALSE)

# split pairs by chromosome
pairs_chrs <- split(pairs, f = pairs[[1]])

# create output directory
dir.create(snakemake@output[[1]], recursive = TRUE)

# save pairs into files per chromosome
for (n in names(pairs_chrs)) {
  write_tsv(pairs_chrs[[n]],
            file = file.path(snakemake@output[[1]], paste0("cre_gene_pairs.", n, ".tsv.gz")),
            col_names = FALSE)
}
