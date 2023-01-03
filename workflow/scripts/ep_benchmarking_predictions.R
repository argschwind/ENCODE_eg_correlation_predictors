# Create E-G predictions from simple correlation coeffcients

save.image("pred.rda")
stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# prepare data -------------------------------------------------------------------------------------

message("Loading data")

# load correlation results
cor <- fread(snakemake@input$cor)

# split cCRE id into coordinates
cor <- separate(cor, col = cre, into = c("chr", "start", "end"), sep = ":|-", remove = FALSE)

# load sCRE and TSS quantifications
cre_quants <- fread(snakemake@input$cres)
tss_quants <- fread(snakemake@input$tss)

# add unique CRE identifier to cre quantifications
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

# create predictions -------------------------------------------------------------------------------

# get quantifications for given cell type
cell_type <- "K562_ID_2644"
cell_type_cre_quants <- cre_quants_mat[, cell_type]
cell_type_cre_quants <- data.table(cre = names(cell_type_cre_quants), reads = cell_type_cre_quants)

# add to correlation data
cor <- merge(cor, cell_type_cre_quants, by = "cre", all.x = TRUE)

# compute correlation * activity predictors
cor <- cor %>% 
  mutate(Score.pearson.activity = pearson * reads,
         Score.spearman.activity = spearman * reads)

# reformat to EPBenchmarking predictions format
output <- cor %>% 
  mutate(class = NA_character_, TargetGeneEnsemblID = NA_character_, TargetGeneTSS = NA_integer_,
         CellType = "K562", DistanceToTSS = NA_integer_) %>% 
  select(chr, start, end, name = cre, class, TargetGene = tss, TargetGeneEnsemblID, TargetGeneTSS,
         CellType, Score.pearson = pearson, Score.spearman = spearman, Score.pearson.activity,
         Score.spearman.activity)

# save to output file
write_tsv(output, file = snakemake@output[[1]])
