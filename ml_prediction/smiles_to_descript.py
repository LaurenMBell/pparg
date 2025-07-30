# 1) clean data into antagonists (IC50) and agonists (EC50)
# 2) convert SMILES code to RDKit descriptors and save those for further use

import pandas as pd
from rdkit import Chem

# read in file
path = "/Users/laurenbell/Desktop/pparg/ml_prediction/bioactivity.tsv"
df = pd.read_csv(path, sep="\t")


# filter to criteria: has values and is what we actually want 
# (EC50 and IC50 > 5) and subset to what we want
df = df[
    (df["Smiles"].notna()) & 
    (df["pChEMBL Value"] >= 5) &
    (df["Standard Type"].isin(["IC50", "EC50"]))
]

df = df[["Smiles", "pChEMBL Value", "Standard Type"]]

df_ag = df[df["Standard Type"] == "EC50"]
df_ant = df[df["Standard Type"] == "IC50"]

# feature selection: 
