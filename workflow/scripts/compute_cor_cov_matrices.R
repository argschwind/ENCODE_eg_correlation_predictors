## Compute correlation and covariance matrices from CRE or TSS quantifications

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

# load quantifications
quant <- fread(snakemake@input[[1]])

# convert quantifications to matrix
quant_mat <- as.matrix(quant[, -1], rownames = quant$id)

# compute correlation among all samples
cor_mat <- cor(quant_mat)
cov_mat <- cov(quant_mat)

# save correation and covariance matrices to files
write.table(cor_mat, file = snakemake@output$cor_mat, sep = "\t")
write.table(cov_mat, file = snakemake@output$cov_mat, sep = "\t")

# remove rownames from matrices for heatmap plots
rownames(cor_mat) <- NULL
rownames(cov_mat) <- NULL

# create heatmaps of correlation and covariance among samples and save to files
pdf(snakemake@output$cor_plot, width = 10, height = 12)
pheatmap(cor_mat)
dev.off()

pdf(snakemake@output$cov_plot, width = 10, height = 12)
pheatmap(cov_mat)
dev.off()
