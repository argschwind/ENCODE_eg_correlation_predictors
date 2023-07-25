## Rules to compute DNase - RNA correlations

# download RNA counts matrix
rule download_rna_counts_table:
  output: "resources/PolyA_RNAseq.TPMs.matrix.gz"
  params:
    url = config["rna_matrix_url"]
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "wget -O {output} {params.url}"

# reformat RNA-seq counts table to required count matrix format
rule reformat_rna_counts_table:
  input: 
    rna = "resources/PolyA_RNAseq.TPMs.matrix.gz",
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz"
  output: "results/RNA/tss_quantifications/RNA_seq/tss_quantifications.allSamples.counts.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/reformat_rna_counts_table.R"