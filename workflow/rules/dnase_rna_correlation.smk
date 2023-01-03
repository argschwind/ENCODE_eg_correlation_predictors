## Rules to compute DNase - RNA correlations

# download gencode v29 annotations
rule download_gencode:
  output: "resources/gencode.v29.annotation.gtf.gz"
  params:
    url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"
  shell:
    "wget -O {output} {params.url}"

# reformat RNA-seq table to required input data
rule reformat_rna_table:
  input: "resource/PolyA_RNAseq.TPMs.matrix.gz"
  output: "results/tss_quantifications/PolyA_RNAseq.counts.tsv.gz"
  conda: 
    
