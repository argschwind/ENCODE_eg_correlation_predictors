## Add samples and read files (bam or TagAlign) from a lookpup table to the snakemake config object 

import pandas as pd

# get file containing samples and read files
samples_file = config['sample_files']['RNA_samples']

# load samples file with sample id as index
samples_table = pd.read_csv(samples_file, sep = '\t', usecols = ['sample', 'DNase_seq', 'RNA_seq'],
                            index_col = 'sample')

# convert to dictionary
samples_dict = samples_table.to_dict('index')

# add samples dictionary to snakemake config object
config['RNA_samples'] = samples_dict
