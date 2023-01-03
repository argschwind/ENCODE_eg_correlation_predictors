## Filter TSS files for TSSs only on regular chromosomes

library(tidyverse)
library(here)

# load input TSS file
tss_file <- here("resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz")
tss <- read_tsv(tss_file, col_names = FALSE, show_col_types = FALSE)

# filter for regular chromosomes
regular_chrs <- paste0("chr", c(seq_len(22), "X", "Y"))
tss <- filter(tss, X1 %in% regular_chrs)

# write output to same file as input
write_tsv(tss, file = tss_file, col_names = FALSE)
