## Rules to compute CRE - TSS correlations using either DNase-only or DNase + RNA-seq

# Create input -------------------------------------------------------------------------------------

# download ENCODE cCREs
rule download_ccres:
  output: temp("resources/GRCh38-cCREs.V4.bed.gz")
  params:
    url = "http://users.wenglab.org/moorej3/Registry-cCREs-WG/V4-Files/GRCh38-cCREs.V4.bed.gz"
  shell:
    "wget -O {output} {params.url}"
    
# sort cCREs according to genomic coordinates and compress
rule sort_ccres:
  input: 
    cres = "resources/GRCh38-cCREs.V4.bed.gz",
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: "resources/GRCh38-cCREs.V4.sorted.bed.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools sort -i {input.cres} -g {input.chrs} | gzip > {output}"
    
# sort TSSs according to genomic coordinates
rule sort_tss:
  input:
    tss  = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz",
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.sorted.bed.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools sort -i {input.tss} -g {input.chrs} | gzip > {output}"
    
# Count DNase-seq reads in cCREs and TSSs ----------------------------------------------------------

# resize cCREs to fixed length (+/- x bp from center)
rule resize_cres:
  input: 
    cres = "resources/GRCh38-cCREs.V4.sorted.bed.gz",
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: "resources/GRCh38-cCREs.V4.sorted.{ext}bp.bed.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "zcat {input.cres} | "
    """awk '{{X={wildcards.ext}; """
    """c=(int($2)+int($3))/2; """
    """printf("%s\\t%d\\t%d\\t%s\\t%s\\t%s\\n",$1,(c-X<0?0:c-X),c+X,$4,$5,$6);}}' | """
    "gzip > {output}"

# count reads in cCREs
rule count_reads_cres:
  input:
    cres = "resources/GRCh38-cCREs.V4.sorted.{ext}bp.bed.gz",
    reads = lambda wildcards: list(str.split(config['abc_samples'][wildcards.sample][wildcards.assay], ",")),
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: temp("results/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.{sample}.counts.bed.gz")
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools coverage -counts -sorted -a {input.cres} -b {input.reads} -g {input.chrs} | "
    "gzip > {output}"
    
# count reads in TSS
rule count_reads_tss:
  input:
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.sorted.bed.gz",
    reads = lambda wildcards: list(str.split(config['abc_samples'][wildcards.sample][wildcards.assay], ",")),
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: temp("results/tss_quantifications/{assay}/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.{sample}.counts.bed.gz")
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools coverage -counts -sorted -a {input.tss} -b {input.reads} -g {input.chrs} | "
    "gzip > {output}"
    
# count reads in cCREs for all ABC samples, combine into one table and add original cCRE coordinates
rule combine_counts_cres:
  input:
    read_quant = expand("results/cre_quantifications/{{assay}}/GRCh38-cCREs.V4.sorted.{{ext}}bp.{sample}.counts.bed.gz",
      sample = config['abc_samples']),
    elements = "resources/GRCh38-cCREs.V4.sorted.bed.gz"
  output: "results/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.allSamples.counts.tsv.gz"
  params:
    type = "cre"
  threads: 5
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/combine_counts.R"
    
# count reads in TSSs for all ABC samples, combine into one table and add original cCRE coordinates
rule combine_counts_tss:
  input:
    read_quant = expand("results/tss_quantifications/{{assay}}/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.{sample}.counts.bed.gz",
      sample = config['abc_samples']),
    elements = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz"
  output: "results/tss_quantifications/{assay}/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.allSamples.counts.tsv.gz"
  params:
    type = "tss"
  threads: 5
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/combine_counts.R"
    
# Compute correlation ------------------------------------------------------------------------------
    
# normalize CRE or TSS quantifications for sequencing depth and log transform
rule normalize_read_counts:
  input: "results/{type}/{assay}/{file}.allSamples.counts.tsv.gz"
  output: "results/{type}/{assay}/{file}.allSamples.counts.normalized.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "12G"
  script:
    "../scripts/normalize_counts.R"

# compute correlation and covariance matrices from quantifications between biosamples
rule compute_cor_cov_matrices:
  input: "results/{type}_quantifications/{assay}/{file}.allSamples.counts.normalized.tsv.gz"
  output: 
    cor_mat = "results/{type}_quantifications/{assay}/{file}.correlation_matrix.tsv",
    cov_mat = "results/{type}_quantifications/{assay}/{file}.covariance_matrix.tsv",
    cor_plot = "results/{type}_quantifications/{assay}/{file}.correlation_plot.pdf",
    cov_plot = "results/{type}_quantifications/{assay}/{file}.covariance_plot.pdf"
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/compute_cor_cov_matrices.R"

# make all enhancer-gene pairs within 1mb
rule make_cre_gene_pairs:
  input:
    cres = "resources/GRCh38-cCREs.V4.sorted.bed.gz",
    tss  = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.sorted.bed.gz"
  output: temp("results/cre_gene_pairs.tsv.gz")
  params:
    max_dist = 1e6
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools window -a {input.cres} -b {input.tss} -w {params.max_dist} | gzip > {output}"

# split pairs into batches
rule split_pairs_batches:
  input: "results/cre_gene_pairs.tsv.gz"
  output:
    temp(expand("results/batches/cre_gene_pairs.batch{batch}.tsv.gz",
      batch = [*range(1, config["batches"] + 1)]))
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/split_pairs_batches.R"

# compute simple correlation
rule compute_correlation:
  input:
    pairs = "results/batches/cre_gene_pairs.batch{batch}.tsv.gz",
    cres = "results/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.allSamples.counts.normalized.tsv.gz",
    tss  = "results/tss_quantifications/{assay}/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.allSamples.counts.normalized.tsv.gz",
    cre_cor = "results/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.correlation_matrix.tsv"
  output: temp("results/{method}/cor_output.{assay}.{ext}bp.batch{batch}.tsv.gz")
  log: "results/logs/compute_correlation.{method}.{assay}.{ext}bp.batch{batch}.log"
  params:
    min_reads = 1
  threads: 24
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "64G",
    time = "12:00:00"
  script:
    "../scripts/compute_correlation.R"
    
# combine results from chromosomes
rule combine_batches:
  input:
    expand("results/{{method}}/cor_output.{{assay}}.{{ext}}bp.batch{batch}.tsv.gz",
      batch = [*range(1, config["batches"] + 1)])
  output: "results/{method}/cor_output.{assay}.{ext}bp.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "16G"
  shell:
    "csvtk concat {input} | gzip > {output}"
    
# Create output ------------------------------------------------------------------------------------ 

# compute all correlation measurements
rule combine_all_correlations:
  input:
    expand("results/{method}/cor_output.{{assay}}.{{ext}}bp.tsv.gz",
      method = ["pearson", "spearman", "gls"])
  output: "results/correlation.{assay}.{ext}bp.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "48G"
  script:
    "../scripts/combine_correlations.R"
    
# # create E-P benchmarking prediction files
# rule create_predictions:
#   input:
#     cres = "results/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.allSamples.counts.bed.gz",
#     tss  = "results/tss_quantifications/{assay}/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.allSamples.counts.bed.gz",
#     cor  = expand("results/{method}/cor_output.{{assay}}.{{ext}}bp.tsv.gz",
#                   method = ['pearson', 'spearman', 'gls'])
#   output: "results/Correlation.predictions.{assay}.{ext}bp.tsv.gz"
#   conda: "../envs/cre_correlation_predictors.yml"
#   resources:
#     mem = "48G"
#   script:
#     "../scripts/ep_benchmarking_predictions.R"
