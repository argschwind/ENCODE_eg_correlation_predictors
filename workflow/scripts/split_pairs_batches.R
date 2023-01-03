## Split cCRE-gene pairs into batches

# required packages
suppressPackageStartupMessages({
  library(data.table)
})

# load cCRE-gene pairs
pairs <- fread(snakemake@input[[1]], header = FALSE)

# number of batches
n_batches <- length(snakemake@output)

# split pairs into more or less equally sized batches
rows_batches <- sort(seq_len(nrow(pairs)) %% n_batches)
batches <- split(pairs, f = rows_batches)

# save pairs into files per chromosome
for (i in seq_along(batches)) {
  fwrite(batches[[i]], file = snakemake@output[[i]], col.names = FALSE, sep = "\t")
}
