## Add metadata for all bam files to config object. Requires a table with metadata for bam file
## accessions

# import modules
import pandas as pd

# get all run type files
meta_files = config['bam_metadata']
meta_files = list(meta_files.values())

# load all run type files and store in list
metadata = []
for file in meta_files:
  metadata_one_file = pd.read_csv(file, sep = '\t')
  metadata.append(metadata_one_file)
  
# combine into one table
metadata = pd.concat(metadata)
metadata = metadata.drop_duplicates()

# make dictionary and add to config object
metadata = metadata.set_index('file')
metadata_dict = metadata.to_dict('index')

# add to config object
config['bam_metadata'] = metadata_dict
