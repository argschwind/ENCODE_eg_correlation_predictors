# ENCODE4 distal regulation E-G correlation predictors
Compute simple enhancer-gene regulatory predictions based on correlation of chromatin accessibility
at cCREs and accessibility at gene promoters or RNA expression levels of genes across different
biosamples.

# Computing DNase-DNase correlations
To compute pearson and spearman correlation coefficients and Generalized Least Square regression
coefficients between DNase-seq signal at cCREs and DNase-seq signal at promoters within 1Mb, run the
following commands. This also installs all dependencies using conda and downloads DNase-seq bam
files from the ENCODE portal. Use the `scratch` parameter in the `config.yml` to specify where
temporary bam files should be stored.

Because of the large number of correlations that are computed and the memory and disk space
requirements, the current implementation of the workflow is only feasible to be executed on an HPC
system.

```sh
# edit 'scratch' parameter to specify location to store bam files, e.g. scratch
vim config/config.yml

# compute all DNase-DNase correlation using parallelization (-n = dryrun, remove to run). adapt this
# command for submission to your cluster based on snakemake documentation
snakemake --use-conda results/DNase/correlation.DNase_seq-DNase_seq.100bp.tsv.gz -j100 -n
```

# Computing DNase-RNA correlations
The workflow can also compute correlations between DNase-seq signal at enhancers with RNA expression
of protein-coding genes with promoters within 1Mb. RNA expression is taken from ENCODE RNA-seq gene
quantification files. To download these files and compute DNase-RNA correlations simply run the
following commands:

```sh
# edit 'scratch' parameter to specify location to store bam files, e.g. scratch
vim config/config.yml

# download DNase-seq and RNA-seq data, and compute DNase-RNA correlations
snakemake --use-conda results/RNA/correlation.DNase_seq-RNA_seq.100bp.tsv.gz -j100 -n
```

# Computing correlation metrics for other sample lists
To compute DNase-DNase or DNase-RNA correlations for other sets of samples, add ENCODE bam and
gene quantification accession IDs to either `config/dnase_dnase_corr_input_files.tsv` or 
`config/dnase_rna_corr_input_files.tsv`. New sample sets can also be provided in new files following
the same format, which can be added to the workflow in the `config/config.yml`. For DNase-seq bam
files the workflow also requires information whether the reads are single- or paired-ended in 
`config/dnase_dnase_corr_dnase_bam_metadata.tsv` or `dnase_rna_corr_dnase_bam_metadata.tsv`,
respectively in a new file following the same format.
