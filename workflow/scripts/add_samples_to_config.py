## Add samples and read files (bam or TagAlign) from a lookpup table to the snakemake config object 

# import modules
import pandas as pd

# function to process one sample list
def process_sample_list(name, file):
  samples_table = pd.read_csv(file, sep = '\t', index_col = 'sample')
  samples_dict = samples_table.to_dict('index')
  config[name] = samples_dict

# get all files containing input files
sample_files = config['sample_files']

# get names of all provided sample lists
names = list(sample_files.keys())

# process all sample files and add contents to snakemake config
for n in names:
  process_sample_list(n, sample_files[n])
