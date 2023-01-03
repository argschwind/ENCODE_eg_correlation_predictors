## Combine original element coordinates with read counts across multiple biosamples

# save.image("combine.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(BiocParallel)
})

# register parallel backend if specified (if more than 1 thread provided)
if (snakemake@threads > 1) {
  register(MulticoreParam(workers = snakemake@threads))
} else {
  register(SerialParam())
}

# cre or tss coordinates
coords <- fread(snakemake@input$elements, header = FALSE)

# files containing read counts for all samples
quant_files <- snakemake@input$read_quant
names(quant_files) <- sub(".+bp\\.(.+)\\.counts.+", "\\1", basename(quant_files))

# function to read quantification file and only retain column with read counts
read_quant_file <- function(file) {
  quant <- fread(file, header = FALSE, nThread = 1)
  return(quant[[ncol(quant)]])
}

# load read counts from all quantification files
read_counts <- bplapply(quant_files, FUN = read_quant_file)

# combine read counts into one table
read_counts <- as.data.table(bind_cols(read_counts))

# create element names depending on whether quantifications are for CREs of TSSs
if (snakemake@params$type == "cre") {
  elements <- paste0(coords[[1]], ":", coords[[2]], "-", coords[[3]])
} else if (snakemake@params$type == "tss") {
  elements <- coords[[4]]
} else {
  stop("Invalid type parameter", call. = FALSE)
}

# add elements to read_counts
read_counts <- cbind(id = elements, read_counts)

# write to output
fwrite(read_counts, file = snakemake@output[[1]], sep = "\t")
