## Compute simple correlation coefficients between CRE and TSS activity

save.image("corrTSS.rda")
stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(BiocParallel)
  library(readr)
})

# register parallel backend if specified (if more than 1 thread provided)
if (snakemake@threads > 1) {
  message("Registering parallel backend with ", snakemake@threads, " cores.")
  register(MulticoreParam(workers = snakemake@threads))
} else {
  message("Registering serial backend.")
  register(SerialParam())
}

# prepare data -------------------------------------------------------------------------------------

message("Loading data")

# load E-G pairs
pairs <- fread(snakemake@input$pairs, header = FALSE)

# load cre and tss read counts
cre_quants <- fread(snakemake@input$cres)
tss_quants <- fread(snakemake@input$tss)

# add cre identifiers to pairs and cre read counts
pairs$cre <- paste0(pairs$V1, ":", pairs$V2, "-", pairs$V3)
cre_quants$cre <- paste0(cre_quants$`#chr`, ":", cre_quants$start, "-", cre_quants$end)

# create matrices containing cre and tss reads
cre_quants_mat <- as.matrix(cre_quants[, -c(1:6)], rownames = "cre")
tss_quants_mat <- as.matrix(tss_quants[, -c(1:3, 5, 6)], rownames = "gene")

# normalize read counts ----------------------------------------------------------------------------

message("Normalizing read counts")

# normalize cre read counts to a total of 1e6 reads per biosample
cre_norm_factors <- colSums(cre_quants_mat) / 1e6
cre_quants_mat <- t(t(cre_quants_mat) / cre_norm_factors)

# normalize tss read counts to a total of 1e6 reads per biosample
tss_norm_factors <- colSums(tss_quants_mat) / 1e6
tss_quants_mat <- t(t(tss_quants_mat) / tss_norm_factors)

# convert to log counts
cre_quants_mat <- log1p(cre_quants_mat)
tss_quants_mat <- log1p(tss_quants_mat)

# filter data based on number of reads -------------------------------------------------------------

# filter pairs for cCREs classified as distal enhancers only and create simplified pairs table
pairs <- pairs[grepl(pairs$V6, pattern = "dELS"), ]
pairs <- data.table(cre = pairs$cre, tss = pairs$V10)

# get any cres or TSS with at least a minimum number of reads
cre_mat_filt <- cre_quants_mat[rowSums(cre_quants_mat) >= snakemake@params$min_reads, ]
tss_mat_filt <- tss_quants_mat[rowSums(tss_quants_mat) >= snakemake@params$min_reads, ]

# filter pairs for cres and TSSs passing the filter
pairs <- pairs[pairs$cre %in% rownames(cre_mat_filt) & pairs$tss %in% rownames(tss_mat_filt), ]

# compute correlations -----------------------------------------------------------------------------

message("Computing correlation")

# helper function to try calculating correlation and capture errors and warnings
try_compute_cor <- function(x, y, method) {
  tryCatch(
    withCallingHandlers({
      cor(x, y, method = method)
    }), warning = function(w) {
      message(w)
      invokeRestart("muffleWarning")
    }, error = function(e) {
      message(e)
      return(NA_real_)
    })
}

# function to compute correlation coefficient(s) for one E-G pair
compute_corr <- function(pair, cre_quants, tss_quants) {
  pearson <- try_compute_cor(cre_quants[pair[[1]], ], tss_quants[pair[[2]], ], method = "pearson")
  spearman <- try_compute_cor(cre_quants[pair[[1]], ], tss_quants[pair[[2]], ], method = "spearman")
  data.table(cre = pair[[1]], tss = pair[[2]], pearson = pearson, spearman = spearman)
}

# split pairs into list for parallel lapply
pairs <- asplit(pairs, MARGIN = 1)

# compute correlation across all pairs
corr_pairs <- bplapply(pairs, FUN = compute_corr, cre_quants = cre_mat_filt,
                       tss_quants = tss_mat_filt)

message("Writing output to files")

# combine into one table
output <- rbindlist(corr_pairs)

# save output to file
write_tsv(output, file = snakemake@output[[1]])

message("Done!")

# close log file connection
sink()
sink(type = "message")
close(log)
