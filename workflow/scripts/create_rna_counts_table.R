## Create RNA gene expression level table based on ENCODE RNA-seq quantification files

# required packages
library(data.table)
library(dplyr)
library(tidyr)

# get files containing RNA quantifications
rna_files <- snakemake@input$rna_tsv
names(rna_files) <- tools::file_path_sans_ext(basename(rna_files))

# load all RNA quantification files
rna <- lapply(rna_files, FUN = fread, select = c("gene_id", "TPM"))

# load genes list file
gene_list <- fread(snakemake@input$gene_list)

# load tss annotations
tss <- fread(snakemake@input$tss, col.names = c("chr", "start", "end", "name", "score", "strand"))

# function to remove any version from Ensembl ids in RNA quantifications and filter for genes on
# genes list for one quantification table from one file
filter_rna <- function(rna_table, genes) {
  filt_rna_table <- rna_table %>% 
    filter(!grepl(gene_id, pattern = "ENSG.+\\..+_PAR_.+$")) %>%  # filter out genes on pseudoautosomal regions
    mutate(gene_id = sub("\\..+", "", gene_id)) %>% 
    filter(gene_id %in% genes$gene_id)
  return(filt_rna_table)
}

# filter all RNA quantification tables for genes on gene list
rna <- lapply(rna, FUN = filter_rna, genes = gene_list)

# combine all filtered RNA quantifications into one table in long format
rna <- rbindlist(rna, idcol = "accession")

# re-calculate TPM after filtering
rna <- rna %>% 
  group_by(accession) %>% 
  mutate(TPM = TPM*1e6 / sum(TPM)) %>% 
  ungroup()

# extract a list of RNA accession ids for each sample from config object for the given sample set
set <- snakemake@wildcards$set
accessions <- snakemake@config[[set]]
accessions <- vapply(accessions, FUN = function(x) x[["RNA_seq"]], FUN.VALUE = character(1))
accessions <- data.table(sample = names(accessions), accession = accessions)

# add sample id and gene symbols to RNA quantification table
rna <- rna %>% 
  left_join(accessions, by = "accession") %>% 
  left_join(select(gene_list, -gene_type), by = "gene_id")

# filter for genes in TSS universe
rna <- filter(rna, gene_name %in% tss$name)

# round TPMs to 5 digits for reproducibility
rna <- mutate(rna, TPM = round(TPM, digits = 5))

# select columns for output and convert to wide format
output <- rna %>% 
  select(-c(accession, gene_id)) %>% 
  rename(id = gene_name) %>% 
  pivot_wider(names_from = sample, values_from = TPM)

# write RNA-seq count matrix to output file
fwrite(output, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
