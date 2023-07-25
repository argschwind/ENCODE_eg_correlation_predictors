# ENCODE4 distal regulation E-G correlation predictors
Compute simple enhancer-gene regulatory predictions based on correlation of chromatin accessibility
at cCREs and accessibility at gene promoters or RNA expression levels of genes across different
biosamples.

# Computing DNase-DNase correlations
To compute pearson and spearman correlation coefficients and Generalized Least Square regression
coefficients between DNase signal at cCREs and DNase signal at promoters within 1Mb, run the
following commands. This also installs all dependencies using conda and dowloads DNase-seq bam files
from the ENCODE portal. Use the `scratch` paramter in the `config.yml` to specify where temporary
bam files should be stored.

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

To compute DNase-DNase correlations on other sets of samples, add ENCODE bam file accessions to the
`config/dnase_dnase_corr_input_files.tsv` file and add information whether the reads in these files
are single- or paired-ended in `config/dnase_dnase_corr_dnase_bam_metadata.tsv`.

# Computing DNase-RNA correlations
For computing correlations between DNase signal at cCREs and RNA expression levels of genes within
1Mb, RNA expression levels were collected manually and stored in a separate file. To download this
file and compute DNase-RNA correlations run the following commands:

```sh
# edit 'scratch' parameter to specify location to store bam files, e.g. scratch
vim config/config.yml

# download and reformat RNA-seq data, and compute DNase-RNA correlations
snakemake --use-conda results/RNA/correlation.DNase_seq-RNA_seq.100bp.tsv.gz -j100 -n

# just download and reformat RNA-seq data without computing correlations
snakemake --use-conda results/RNA/tss_quantifications/RNA_seq/tss_quantifications.allSamples.counts.tsv.gz
```

To compute DNase-RNA correlations on additional samples, add DNase-seq bam files to config files
like previously described and create a table of RNA expression levels using the example files above.
