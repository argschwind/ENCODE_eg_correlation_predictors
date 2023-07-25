## Compute correlations between CRE and TSS activity or fit GLS models

# save.image("compute_correlation.rda")
# stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nlme)
  library(BiocParallel)
})

# register parallel backend if specified (if more than 1 thread provided)
if (snakemake@threads > 1) {
  message("\nRegistering parallel backend with ", snakemake@threads, " cores:")
  register(MulticoreParam(workers = snakemake@threads))
  bpparam()
} else {
  message("\nRegistering serial backend:")
  register(SerialParam())
  bpparam()
}

# prepare data -------------------------------------------------------------------------------------

message("\nLoading data")

# load E-G pairs
pairs <- fread(snakemake@input$pairs)

# load normalized CRE and TSS read counts
cre_quants <- fread(snakemake@input$cres)
tss_quants <- fread(snakemake@input$tss) 

# convert CRE and TSS read counts to matrices
cre_quants_mat <- as.matrix(cre_quants[, -1], rownames = cre_quants$id)
tss_quants_mat <- as.matrix(tss_quants[, -1], rownames = tss_quants$id)

# make sure that samples (columns) are the same (and same order) for CRE and TSS quantifications
samples <- intersect(colnames(cre_quants_mat), colnames(tss_quants_mat))
cre_quants_mat <- cre_quants_mat[, samples]
tss_quants_mat <- tss_quants_mat[, samples]

# get any CREs or TSS with at least a minimum number of reads
cre_quants_mat <- cre_quants_mat[rowSums(cre_quants_mat) >= snakemake@params$min_reads, ]
tss_quants_mat <- tss_quants_mat[rowSums(tss_quants_mat) >= snakemake@params$min_reads, ]

# create simplified table of pairs with only CRE and TSS ids
pairs <- data.table(cre = paste0(pairs$V1, ":", pairs$V2, "-", pairs$V3), tss = pairs$V10)

# filter pairs for CREs and TSSs in quantification matrices and passing the minimum reads filter
pairs <- pairs[pairs$cre %in% rownames(cre_quants_mat) & pairs$tss %in% rownames(tss_quants_mat), ]

# load and process correlation matrix if GLS is used as method -------------------------------------

# get method from wildcards
method <- snakemake@wildcards$method

# if gls is specified load cre sample correlation matrix and create correlation structure for gls
if (method == "gls") {
  
  message("Setting up correlation structure for GLS")
  
  # load cre correlation matrix
  cre_cor_mat <- as.matrix(read.table(snakemake@input$cre_cor, header = TRUE, sep = "\t"))
  
  # make sure that correlation matrix is a positive-definite matrix
  cre_cor_mat <- nearPD(cre_cor_mat)$mat
  
  # create correlation structure from sample correlation matrix for GLS
  cor_structure <- corSymm(cre_cor_mat[lower.tri(cre_cor_mat)], fixed = TRUE)
  
}

# define functions to compute correlations or fit regression models --------------------------------

# helper function to try calculating correlation and capture errors and warnings
try_compute_cor <- function(cor_function, pair, cre_quants, tss_quants, ...) {
  tryCatch(
    withCallingHandlers({
      cor_function(pair, cre_quants = cre_quants, tss_quants = tss_quants, ...)
    }), warning = function(w) {
      message(w)
      invokeRestart("muffleWarning")
    }, error = function(e) {
      message(e)
      return(NA_real_)
    })
}

# function to compute pearson or spearman correlation coefficient for one E-G pair
compute_simple_cor <- function(pair, cre_quants, tss_quants, method) {
  cor_coeff <- cor(cre_quants[pair[[1]], ], tss_quants[pair[[2]], ], method = method)
  setNames(data.frame(pair[[1]], pair[[2]], cor_coeff), c("cre", "tss", method))
}

# function to fit GLS for one enhancer-gene pair
compute_gls <- function(pair, cre_quants, tss_quants, cor_structure = NULL) {
  
  # create table with tss and cre quantifications
  quants <- data.frame(tss = tss_quants[pair[[2]], ], cre = cre_quants[pair[[1]], ])
  
  # fit model
  fit <- gls(model = tss ~ cre, data = quants, correlation = cor_structure)
  
  # extract coefficient and p-values to create output
  t_table <- summary(fit)$tTable
  output <- data.table(cre = pair[[1]], tss = pair[[2]], glsCoefficient = t_table["cre", "Value"],
                       glsPvalue  = t_table["cre", "p-value"])
  
  return(output)
  
}

# compute correlation or fit models ----------------------------------------------------------------

# split pairs into list for parallel lapply
pairs <- asplit(pairs, MARGIN = 1)

# compute correlation or fit across all pairs
message("Computing correlations")
if (method %in% c("pearson", "spearman")) {
  bpparam()
  output <- bplapply(pairs, FUN = try_compute_cor, cor_function = compute_simple_cor,
                     cre_quants = cre_quants_mat, tss_quants = tss_quants_mat, method = method)
  
} else if (method == "gls") {
  
  output <- bplapply(pairs, FUN = try_compute_cor, cor_function = compute_gls,
                     cre_quants = cre_quants_mat, tss_quants = tss_quants_mat,
                     cor_structure = cor_structure)
  
} else {
  
  stop("Invalid method wildard value: ", method,
       ". Needs to be one of 'pearson', 'spearman' or 'gls'")
  
}

# combine output list into one table
output <- rbindlist(output)

# save output to file
message("Writing results to output file")
fwrite(output, file = snakemake@output[[1]], sep = "\t")

message("Done!")

# close log file connection
sink()
sink(type = "message")
close(log)
