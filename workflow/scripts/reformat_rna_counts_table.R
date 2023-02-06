## Reformat RNA-seq counts table to required count matrix format

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load input RNA-seq counts table
rna <- fread(snakemake@input$rna)

# load tss annotations
tss <- fread(snakemake@input$tss, col.names = c("chr", "start", "end", "name", "score", "strand"))

# reformat for output
rna <- rna %>% 
  select(-c(gene_id, chrom, TSS)) %>% 
  rename(id = symbol)

# only retain data on genes also found in used TSS annotations
rna <- filter(rna, id %in% unique(tss$name))

# write reformatted count matrix to output file
fwrite(rna, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
