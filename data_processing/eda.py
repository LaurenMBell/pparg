"""
4/16/26

Tutorials/books referenced: 
 - https://medium.com/@post.gourang/a-holistic-guide-to-exploratory-data-analysis-eda-for-machine-learning-and-deep-learning-bc4f18f0143b
 - https://www.linkedin.com/pulse/best-python-libraries-exploratory-data-analysis-eda-suraj-kumar-soni-pqruc/
 - http://ciml.info/dl/v0_99/ciml-v0_99-ch05.pdf
  - https://darkomedin-datascience.medium.com/data-science-for-drug-discovery-research-morgan-fingerprints-using-alanine-and-testosterone-92a2c69dd765
 
What is the problem I'm aiming to solve here? 
    I want to classify molecule-molecule interactions using random forests. Or, 
    I guess put more simply, I want to classify bioactivity (IC50) of untested
    compounds using the activity data of tested compounds.

Features for model (standard type and SMILES code don't count): 
1) Molecular weight: continuous 
2) #RO5 volations: categorical, 4 categories, 0-3 range
3) HA (# of heavy atoms): continuous, 
4) HBD (h-bond donors): categorical
5) HBA (h-bond acceptors): categ7orical
6) LogP: continuous
7) Morgan fingerprints: continuous
"""

import matplotlib.pyplot as plt
import pandas as pd 
import sweetviz as sv
from featSelection import Features

def summary_stats(df):    
    report = sv.analyze([df, 'Train'], target_feat = "pChEMBL Value")
    report.show_html('reports/alldata_report.html')

def feature_stats(df, label):
    mfps = [col for col in df.columns if str(col).startswith("mfp_")]
    df_clean = df.drop(columns=mfps) #drop all the mfp features and smiles codes
    #df_clean = df_clean.drop(["Smiles"])
    report = sv.analyze([df_clean, 'Train'], target_feat = "activity")
    report.show_html(f"reports/{label}_report.html")


def main():
    bioactivity = pd.read_csv("data/bioactivity.tsv", delimiter="\t")
    bioactivity = bioactivity.dropna(subset="pChEMBL Value") #drop any NaN IC50/EC50 values

    summary_stats(bioactivity) # create EDA report with sweetviz

    # filter to criteria
    cleaned_bioactivity = bioactivity[ 
        (bioactivity["Smiles"].notna()) & # need for RDKit 
        (bioactivity["pChEMBL Value"].notna()) & # target feature
        (bioactivity["Standard Type"].isin(["IC50", "EC50"]))] #only activators and inhibitors

    feature_data = cleaned_bioactivity[["Smiles", "pChEMBL Value", "Standard Type"]]
    feature_data['activity'] = feature_data['pChEMBL Value'].apply(Features.activity) #map to binary activity labels

    #split into agonist and antagonist datasets
    df_ag = feature_data[feature_data["Standard Type"] == "EC50"]
    df_ant = feature_data[feature_data["Standard Type"] == "IC50"]

    # feature selection
    ag_feats = df_ag["Smiles"].apply(lambda x: pd.Series(Features.calc_descriptors(x)))
    df_ag = pd.concat([df_ag, ag_feats], axis = 1)
    ant_feats = df_ant["Smiles"].apply(lambda x: pd.Series(Features.calc_descriptors(x)))
    df_ant = pd.concat([df_ant, ant_feats], axis = 1)

    df_ag.to_csv("data/agonists.csv", index=False)
    df_ant.to_csv("data/antagonists.csv", index=False)

    #get EDA for ag and ant feature datasets
    feature_stats(df_ag, "agonists")
    feature_stats(df_ant, "antagonists")

if __name__ == "__main__":
    main()