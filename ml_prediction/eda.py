# Exploratory data analysis of the RDKit descriptors!!

import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd

ag = pd.read_csv("ml_prediction/data/agonists.csv")
ant = pd.read_csv("ml_prediction/data/antagonists.csv")

# show distribution of molecular weight
def mw_dist(data, title):
    fig, ax = plt.subplots()

    ax.hist(data["molecular_weight"], bins=50)
    plt.axvline(x=500, color="red", label="Lipinksi Threshold")
    plt.legend()

    plt.xlabel("Average Molecular Weight")
    plt.ylabel("Frequency")
    plt.title(title)

    plt.show()

#mw_dist(ag, "Distribution of Molecular Weight in PPARG Agonists")
#mw_dist(ant,"Distribution of Molecular Weight in PPARG Antagonists")

# show # of HBA and HBDs in both ag and ant data 
# (I've learned this is bad visualization and tried histograms instead)
def hb_scatter(data1, data2, title):
    fig, ax = plt.subplots()

    ax.scatter(data1["HBA"], data1["HBD"], color="red")
    ax.scatter(data2["HBA"], data2["HBD"], color="blue")

    plt.axvline(x=10, color="red", label="Lipinksi Threshold for HBA")
    plt.axhline(y=5, color="orange", label="Lipinksi Threshold for HBD")


    plt.legend()
    plt.xlabel("Number of H-bond Acceptors")
    plt.ylabel("Number of H-bond Donors")
    plt.title(title)

    plt.show()

#hba_scatter(ag, ant, "H-Bond Acceptors/Donors Across PPARG\nAntagonists and Agonists")

# distribution of hba stacked for both datasets
def hba_hist(data1, data2, title):
    fig, ax = plt.subplots()
    ax.hist([data1["HBA"], data2["HBA"]], bins = 10, stacked=True, color = ["blue", "green"])
    
    plt.axvline(x=10, color="red", label="Lipinksi Threshold for HBA")

    plt.legend()
    plt.xlabel("Number of H-bond Acceptors")
    plt.ylabel("Frequency")
    plt.title(title)

    plt.show()

#hba_hist(ag, ant, Distribution of H-Bond Acceptors Across\nPPARG Antagonists and Agonists"))


# distribution of hbd stacked for both datasets
def hbd_hist(data1, data2, title):
    fig, ax = plt.subplots()
    ax.hist([data1["HBD"], data2["HBD"]], bins = 10, stacked=True, color = ["blue", "orange"])
    
    plt.axvline(x=5, color="red", label="Lipinksi Threshold for HBD")

    plt.legend()
    plt.xlabel("Number of H-bond Donors")
    plt.ylabel("Frequency")
    plt.title(title)

    plt.show()

#hbd_hist(ag, ant, "Distribution of H-Bond Donors Across\nPPARG Antagonists and Agonists")

