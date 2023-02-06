## Rules to compute DNase - RNA correlations

ruleorder: sort_bam > download_bam

# reformat RNA-seq counts table to required count matrix format
rule reformat_rna_counts_table:
  input: 
    rna = "resources/PolyA_RNAseq.TPMs.matrix.gz",
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz"
  output: "results/RNA_samples/tss_quantifications/RNA_seq/tss_quantifications.allSamples.counts.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/reformat_rna_counts_table.R"
  
# download bam file from ENCODE portal
rule download_bam:
  output: temp(config["scratch"] + "/correlation/bam/{accession}.bam")
  params:
    base_url = "https://www.encodeproject.org/files"
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "wget -O {output} {params.base_url}/{wildcards.accession}/@@download/{wildcards.accession}.bam"
    
# sort and index bam files
rule sort_bam:
  input: config["scratch"] + "/correlation/bam/{accession}.bam"
  output:
    bam = config["scratch"] + "/correlation/bam/{accession}.sorted.bam",
    bai = config["scratch"] + "/correlation/bam/{accession}.sorted.bam.bai"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "32G",
    time = "4:00:00"
  shell:
    "samtools sort -o {output.bam} {input}; samtools index {output.bam}"
