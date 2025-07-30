# Random Forest model to predict the IC50 score for antagonists!

import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen

import argparse

# same SMILES to descriptors conversion
def calc_descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    
    # lipinski features
    mw = Descriptors.MolWt(m)
    hba = Lipinski.NumHAcceptors(m)
    hbd = Lipinski.NumHDonors(m)
    logp = Crippen.MolLogP(m)

    return mw, hba, hbd, logp

# prediction logic!! input SMILES code, get activity prediction
def predict_from_smiles(smiles, model):
    descriptors = calc_descriptors(smiles)
    ddf = pd.DataFrame()
    ddf = pd.DataFrame([descriptors], 
                    columns=['molecular_weight', 'HBA', 'HBD', 'LogP'])
    
    prediction = model.predict(ddf)[0]
    prob = model.predict_proba(ddf)[0]

    if prediction == 1:
        return 'active', max(prob), descriptors
    else:
        return 'inactive', max(prob), descriptors

#===========================================================================

def main():
    #set up input args
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles")

    args = parser.parse_args()

    # reading in and mapping 
    df = pd.read_csv("data/antagonists.csv")
    df['activity'] = df['activity'].map({'inactive' : 0, 'active' : 1})

    # splitting the data
    X = df[['molecular_weight', 'HBA', "HBD", "LogP"]]
    y = df['activity']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

    #fit model with training data
    rf = RandomForestClassifier()
    rf.fit(X_train, y_train)

    #testing the model 
    y_pred = rf.predict(X_test)

    # evaluating accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy: ", accuracy)



    result = predict_from_smiles(args.smiles, rf)
    print(f"Input code: {args.smiles}")
    print(f"Prediction: {result}")

if __name__ == "__main__":
    main()



