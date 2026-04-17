"""
1) clean data into antagonists (IC50) and agonists (EC50)
2) convert SMILES code to RDKit descriptors and save those for further use

Made it a class so I could use it in eda.py

Tutorials referenced: 
 - https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html#:~:text=I%20already%20showed%20how%20to,(339%2C%20984)
 - https://darkomedin-datascience.medium.com/data-science-for-drug-discovery-research-morgan-fingerprints-using-alanine-and-testosterone-92a2c69dd765

IC50/EC50 threshold choice: 
https://www.researchgate.net/post/Acceptable-IC50-drug-concentration-for-MTT-essay#:~:text=In%20most%20cases%2C%20the%20IC50,and%20excretion%20in%20our%20body.

The value provided by ChEMBL is -log(molar IC50), so I just 
converted it from there to 5.5. pChEMBL >= 5/5 is active, and 
pChEMBL < 5.5 is inactive. This is admittedly a pretty weak
threshold, but for a project like this its fine. I wanted a 
balanced anough dataset for the models.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdFingerprintGenerator
import numpy as np

class Features:
    @staticmethod

    # function for getting from smiles code -> features we want
    def calc_descriptors(smiles):
        m = Chem.MolFromSmiles(smiles)
        
        # lipinski features
        mw = Descriptors.MolWt(m)
        hba = Lipinski.NumHAcceptors(m)
        hbd = Lipinski.NumHDonors(m)
        logp = Crippen.MolLogP(m)

        bit={}
        #theoretically this is depreciated but I don't really care rn
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        fp = mfpgen.GetFingerprint(m)

        ha = m.GetNumHeavyAtoms()

        return mw, hba, hbd, logp, fp, ha

    def activity(val):
        if val >= 5.5:
            return True
        else:
            return False

