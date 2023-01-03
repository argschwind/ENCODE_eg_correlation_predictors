## Normalize CRE or TSS quantifications for sequencing depth and log transform

# required packages
suppressPackageStartupMessages({
  library(data.table)
})

# load raw read counts
counts <- fread(snakemake@input[[1]])

# normalize read counts to a total of 1e6 reads per biosample and log transform
norm_factors <- colSums(counts[, -1]) / 1e6
norm_counts <- sweep(counts[, -1], MARGIN = 2, STATS = norm_factors, FUN = "/")
norm_counts <- log1p(norm_counts)
norm_counts <- cbind(counts[, 1], norm_counts)  # find better solution...

# write normalized data to output file
fwrite(norm_counts, file = snakemake@output[[1]], sep = "\t")
