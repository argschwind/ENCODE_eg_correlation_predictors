## Rules to compute DNase - RNA correlations

# make sure that count table for RNA is produced from combining tsv files
ruleorder: create_rna_counts_table > combine_counts_tss

# input function to get all RNA quantification files for a given sample set
def get_rna_tsv_files(wildcards):
  samples = config[wildcards.set]
  accessions = [value.get('RNA_seq', None) for value in samples.values()]
  file_pattern = config['scratch'] + '/correlation/rna_tsv/{}.tsv'
  files = list(map(file_pattern.format, accessions))
  return(files)

# download RNA quantification files (.tsv)
rule download_rna_tsv:
  output: config["scratch"] + "/correlation/rna_tsv/{accession}.tsv"
  params:
    base_url = "https://www.encodeproject.org/files"
  conda: "../envs/cre_correlation_predictors.yml"
  shell:
    "wget -O {output} {params.base_url}/{wildcards.accession}/@@download/{wildcards.accession}.tsv"
    
# create RNA counts matrix
rule create_rna_counts_table:
  input: 
    rna_tsv = get_rna_tsv_files,
    gene_list = "resources/rnaseq_protein_coding_genes.tsv.gz",
    tss = "resources/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed.gz"
  output: "results/{set}/tss_quantifications/RNA_seq/tss_quantifications.allSamples.counts.tsv.gz"
  conda: "../envs/cre_correlation_predictors.yml"
  script:
    "../scripts/create_rna_counts_table.R"
