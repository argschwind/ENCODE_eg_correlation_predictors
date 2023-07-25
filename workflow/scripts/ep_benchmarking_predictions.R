## Create E-G prediction files from E-G correlation metrics

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

# load cCRE annotations
cres <- fread(snakemake@input$cres, drop = 4:5, col.names = c("chr", "start", "end", "class"))

message("Creating prediction files")

# create unique cCRE identifiers
cres <- mutate(cres, name = paste0(chr, ":", start, "-", end))

# create enhancer id
cor <- mutate(cor, name = paste0(chr, ":", start, "-", end))

# add TSS coordinates and cCRE class to correlation results
cor <- cor %>% 
  left_join(select(tss, TargetGene = V4, TargetGeneTSS = V2), by = "TargetGene") %>% 
  left_join(select(cres, name, class), by = "name")

# compute distance to tss
cor <- cor %>% 
  mutate(enh_center = (end + start) / 2) %>% 
  mutate(DistanceToTSS = round(abs(TargetGeneTSS - enh_center)))

# reformat to EPBenchmarking predictions format including renaming score columns
output <- cor %>% 
  mutate(TargetGeneEnsemblID = NA_character_, CellType = NA_character_) %>% 
  select(chr, start, end, name, class, TargetGene, TargetGeneEnsemblID, TargetGeneTSS,
         CellType, any_of(c(pearson.Score = "pearson", spearman.Score = "spearman",
         glsCoefficient.Score = "glsCoefficient")), DistanceToTSS)

# save to output file
write_tsv(output, file = snakemake@output[[1]])

message("Done")
