# 1) clean data
# 2) convert SMILES code to RDKit descriptors 
#    and save those for further use

import pandas as pd

path = "/pparg_bioactivity.tsv"

df = pd.read_csv(path, sep="\t")



