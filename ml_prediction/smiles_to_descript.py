# 1) clean data into antagonists (IC50) and agonists (EC50)
# 2) convert SMILES code to RDKit descriptors and save those for further use

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen

# function for getting from smiles code -> features we want
def calc_descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    
    # lipinski features
    mw = Descriptors.MolWt(m)
    hba = Lipinski.NumHAcceptors(m)
    hbd = Lipinski.NumHDonors(m)
    logp = Crippen.MolLogP(m)

    return mw, hba, hbd, logp

# read in file
path = "ml_prediction/data/bioactivity.tsv"
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

# feature selection
descript = df_ag["Smiles"].apply(calc_descriptors)
df_ag[["molecular_weight", "HBA", "HBD", "LogP"]] = pd.DataFrame(descript.tolist(), index=df_ag.index)

descript = df_ant["Smiles"].apply(calc_descriptors)
df_ant[["molecular_weight", "HBA", "HBD", "LogP"]] = pd.DataFrame(descript.tolist(), index=df_ant.index)

df_ag.to_csv("ml_prediction/data/agonists.csv", index=False)
df_ant.to_csv("ml_prediction/data/antagonists.csv", index=False)
