# Create E-G prediction files from simple correlation coeffcients

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
})

message("Loading data")

# load correlation results
cor <- fread(snakemake@input$cor)

# load TSS annotations
tss <- fread(snakemake@input$tss, header = FALSE)

# add TSS coordinates to correlation results
cor <- tss %>% 
  select(TargetGene = V4, TargetGeneTSS = V2) %>% 
  left_join(cor, ., by = "TargetGene")

# create enhancer id
cor <- mutate(cor, name = paste0(chr, ":", start, "-", end))

# compute distance to tss
cor <- cor %>% 
  mutate(enh_center = (end + start) / 2) %>% 
  mutate(DistanceToTSS = round(abs(TargetGeneTSS - enh_center)))

# reformat to EPBenchmarking predictions format
output <- cor %>% 
  mutate(class = NA_character_, TargetGeneEnsemblID = NA_character_, CellType = NA_character_) %>% 
  select(chr, start, end, name, class, TargetGene, TargetGeneEnsemblID, TargetGeneTSS,
         CellType, pearson.Score = pearson, spearman.Score = spearman,
         glsCoefficient.Score = GLScoefficient, DistanceToTSS)

# save to output file
write_tsv(output, file = snakemake@output[[1]])
