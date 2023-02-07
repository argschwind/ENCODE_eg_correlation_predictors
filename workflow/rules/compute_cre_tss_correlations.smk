## Rules to compute CRE - TSS correlations using either DNase-only or DNase + RNA-seq

# Python input functions ---------------------------------------------------------------------------

# get mapped reads file (e.g. tagAlign or bam)
def get_read_files(wildcards):
  type = config["read_files"][wildcards.set]
  files = list(str.split(config[wildcards.set][wildcards.sample][wildcards.assay], ","))
  if type == "accession":
    file_pattern = config['scratch'] + '/correlation/bam/{}.sorted.bam'
    files = list(map(file_pattern.format, files))
  return files
  
# Download bam files from ENCODE portal ------------------------------------------------------------

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
    reads = get_read_files,
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: temp("results/{set}/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.{sample}.counts.bed.gz")
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools coverage -counts -sorted -a {input.cres} -b {input.reads} -g {input.chrs} | "
    "gzip > {output}"
    
# count reads in TSS
rule count_reads_tss:
  input:
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.sorted.bed.gz",
    reads = get_read_files,
    chrs = "resources/GRCh38_EBV.chrom.sizes.tsv"
  output: temp("results/{set}/tss_quantifications/{assay}/tss_quantifications.{sample}.counts.bed.gz")
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools coverage -counts -sorted -a {input.tss} -b {input.reads} -g {input.chrs} | "
    "gzip > {output}"
    
# count reads in cCREs for all ABC samples, combine into one table and add original cCRE coordinates
rule combine_counts_cres:
  input:
    read_quant = lambda wildcards:
      expand("results/{{set}}/cre_quantifications/{{assay}}/GRCh38-cCREs.V4.sorted.{{ext}}bp.{sample}.counts.bed.gz",
        sample = config[wildcards.set]),
    elements = "resources/GRCh38-cCREs.V4.sorted.bed.gz"
  output: "results/{set}/cre_quantifications/{assay}/GRCh38-cCREs.V4.sorted.{ext}bp.allSamples.counts.tsv.gz"
  params:
    type = "cre"
  threads: 5
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/combine_counts.R"

# count reads in TSSs for all ABC samples, combine into one table and add original cCRE coordinates
rule combine_counts_tss:
  input:
    read_quant = lambda wildcards:
      expand("results/{{set}}/tss_quantifications/{{assay}}/tss_quantifications.{sample}.counts.bed.gz",
        sample = config[wildcards.set]),
    elements = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz"
  output: "results/{set}/tss_quantifications/{assay}/tss_quantifications.allSamples.counts.tsv.gz"
  params:
    type = "tss"
  threads: 5
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/combine_counts.R"

# Compute correlation ------------------------------------------------------------------------------
    
# normalize CRE or TSS quantifications for sequencing depth and log transform
rule normalize_read_counts:
  input: "results/{set}/{type}/{assay}/{file}.allSamples.counts.tsv.gz"
  output: "results/{set}/{type}/{assay}/{file}.allSamples.counts.normalized.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "12G"
  script:
    "../scripts/normalize_counts.R"

# compute correlation and covariance matrices from quantifications between biosamples
rule compute_cor_cov_matrices:
  input: "results/{set}/{type}_quantifications/{assay}/{file}.allSamples.counts.normalized.tsv.gz"
  output: 
    cor_mat = "results/{set}/{type}_quantifications/{assay}/{file}.correlation_matrix.tsv",
    cov_mat = "results/{set}/{type}_quantifications/{assay}/{file}.covariance_matrix.tsv",
    cor_plot = "results/{set}/{type}_quantifications/{assay}/{file}.correlation_plot.pdf",
    cov_plot = "results/{set}/{type}_quantifications/{assay}/{file}.covariance_plot.pdf"
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/compute_cor_cov_matrices.R"

# make all enhancer-gene pairs within 1mb
rule make_cre_gene_pairs:
  input:
    cres = "resources/GRCh38-cCREs.V4.sorted.bed.gz",
    tss  = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.sorted.bed.gz"
  output: temp("results/cre_gene_pairs/cre_gene_pairs.tsv.gz")
  params:
    max_dist = 1e6
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "bedtools window -a {input.cres} -b {input.tss} -w {params.max_dist} | gzip > {output}"

# split pairs into batches
rule split_pairs_batches:
  input: "results/cre_gene_pairs/cre_gene_pairs.tsv.gz"
  output:
    temp(expand("results/cre_gene_pairs/cre_gene_pairs.batch{batch}.tsv.gz",
      batch = [*range(1, config["batches"] + 1)]))
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/split_pairs_batches.R"

# compute simple correlation
rule compute_correlation:
  input:
    pairs = "results/cre_gene_pairs/cre_gene_pairs.batch{batch}.tsv.gz",
    cres = "results/{set}/cre_quantifications/{cre_assay}/GRCh38-cCREs.V4.sorted.{ext}bp.allSamples.counts.normalized.tsv.gz",
    tss  = "results/{set}/tss_quantifications/{tss_assay}/tss_quantifications.allSamples.counts.normalized.tsv.gz",
    cre_cor = "results/{set}/cre_quantifications/{cre_assay}/GRCh38-cCREs.V4.sorted.{ext}bp.correlation_matrix.tsv"
  output: temp("results/{set}/{method}/cor_output.{cre_assay}-{tss_assay}.{ext}bp.batch{batch}.tsv.gz")
  log: "results/{set}/logs/compute_correlation.{method}.{cre_assay}-{tss_assay}.{ext}bp.batch{batch}.log"
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
    expand("results/{{set}}/{{method}}/cor_output.{{cre_assay}}-{{tss_assay}}.{{ext}}bp.batch{batch}.tsv.gz",
      batch = [*range(1, config["batches"] + 1)])
  output: "results/{set}/{method}/cor_output.{cre_assay}-{tss_assay}.{ext}bp.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "16G"
  shell:
    "csvtk concat {input} | gzip > {output}"
    
# Create output ------------------------------------------------------------------------------------ 

# compute all correlation measurements
rule combine_all_correlations:
  input:
    expand("results/{{set}}/{method}/cor_output.{{cre_assay}}-{{tss_assay}}.{{ext}}bp.tsv.gz",
      method = ["pearson", "spearman", "gls"])
  output: "results/{set}/correlation.{cre_assay}-{tss_assay}.{ext}bp.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "48G"
  script:
    "../scripts/combine_correlations.R"
    
# create E-P benchmarking prediction files
rule create_prediction_files:
  input:
    cres = "resources/GRCh38-cCREs.V4.sorted.bed.gz",
    tss  = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.sorted.bed.gz",
    cor = "results/{set}/correlation.{cre_assay}-{tss_assay}.{ext}bp.tsv.gz"
  output: "results/{set}/correlation.{cre_assay}-{tss_assay}.{ext}bp.predictions.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  resources:
    mem = "12G"
  script:
    "../scripts/ep_benchmarking_predictions.R"
