# Exploratory data analysis of the RDKit descriptors!!

import matplotlib.pyplot as plt
import pandas as pd

red = '#BF0603'
blue = '#26A3D1'
dark_blue = "#124E78"
green = "#679436"
yellow = "#F7CB15"
 


ag = pd.read_csv("ml_prediction/data/agonists.csv")
ag.name = "ag"
ant = pd.read_csv("ml_prediction/data/antagonists.csv")
ant.name = "ant"

# show distribution of molecular weight
def mw_dist(data, title):
    fig, ax = plt.subplots()
    ax.grid(zorder=0)
    ax.hist(data["molecular_weight"], bins=75, color=blue, zorder=3)
    plt.axvline(x=500, color=red, label="Lipinksi Threshold", zorder=4)
    plt.legend()

    plt.xlabel("Average Molecular Weight")
    plt.ylabel("Frequency")
    plt.title(title)

    plt.savefig("ml_prediction/eda_figs/" + data.name + "_hba_dist.png", dpi=300)

mw_dist(ag, "Distribution of Molecular Weight in PPARG Agonists")
mw_dist(ant,"Distribution of Molecular Weight in PPARG Antagonists")

# show # of HBA and HBDs in both ag and ant data 
# (I've learned this is bad visualization and tried histograms instead)
def hb_scatter(data1, data2, title):
    fig, ax = plt.subplots()

    ax.scatter(data1["HBA"], data1["HBD"], color=green)
    ax.scatter(data2["HBA"], data2["HBD"], color=blue)

    plt.axvline(x=10, color=red, label="Lipinksi Threshold for HBA")
    plt.axhline(y=5, color=yellow, label="Lipinksi Threshold for HBD")


    plt.legend()
    plt.xlabel("Number of H-bond Acceptors")
    plt.ylabel("Number of H-bond Donors")
    plt.title(title)

    plt.savefig("ml_prediction/eda_figs/hb_scatter.png", dpi=300)

hb_scatter(ag, ant, "H-Bond Acceptors/Donors Across PPARG\nAntagonists and Agonists")

# distribution of hba side by side for both datasets
def hba_hist(data1, data2, title):
    fig, ax = plt.subplots()
    ax.hist([data1["HBA"], data2["HBA"]], bins = 10, color = [blue, yellow], label=["Agonists", "Antagonists"])
    
    plt.axvline(x=10, color="red", label="Lipinksi Threshold")

    plt.legend()
    plt.xlabel("Number of H-bond Acceptors")
    plt.ylabel("Frequency")
    plt.title(title)

    plt.savefig("ml_prediction/eda_figs/hba_dist.png", dpi=300)

hba_hist(ag, ant, "Distribution of H-Bond Acceptors Across\nPPARG Antagonists and Agonists")


# distribution of hbd side by side for both datasets
def hbd_hist(data1, data2, title):
    fig, ax = plt.subplots()
    ax.hist([data1["HBD"], data2["HBD"]], bins = 10, color = [blue, yellow], label=["Agonists", "Antagonists"])
    
    plt.axvline(x=5, color="red", label="Lipinksi Threshold")

    plt.legend()
    plt.xlabel("Number of H-bond Donors")
    plt.ylabel("Frequency")
    plt.title(title)

    plt.savefig("ml_prediction/eda_figs/hbd_dist.png", dpi=300)

hbd_hist(ag, ant, "Distribution of H-Bond Donors Across\nPPARG Antagonists and Agonists")

