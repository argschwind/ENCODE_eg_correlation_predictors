## Combine correlation results from different methods

# required packages
library(data.table)
library(dplyr)
library(tidyr)

# files containing cre-promoter correlation results
result_files <- unlist(snakemake@input)
names(result_files) <- basename(dirname(result_files))

# read all result files
results <- lapply(result_files, FUN = fread, header = TRUE)

# combine into one data frame
results <- Reduce(function(df1, df2) merge(df1, df2, by = c("cre", "tss"), all = TRUE),
                  results)

# split correlation cre id into coordinates
results <- separate(results, cre, into = c("chr", "start", "end"), sep = ":|-", convert = TRUE)

# select and rearrange columns for output
results <- results %>% 
  select(chr, start, end, TargetGene = tss,
         any_of(c("pearson", "spearman", "glsCoefficient", "glsPvalue")))

# save to output file
fwrite(results, file = snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
