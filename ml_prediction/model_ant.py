# Random Forest model to predict the IC50 score for antagonists!

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from rdkit import Chem
from scipy.stats import randint
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
    # set up input args for the bash script
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

    param_dist = {'n_estimators': randint(50,500), 
                  'max_depth': randint(1,20)}

    # process of optimizing hyperparameters! I found they were:
    # {'max_depth': 6, 'n_estimators': 304}
    
    #rf = RandomForestClassifier()
    #rand_search = RandomizedSearchCV(rf, 
                                     #param_distributions= param_dist,
                                     #n_iter=5,
                                     #cv=5)
    #rand_search.fit(X_train, y_train)

    # use the best of the five models created
    #best_rf = rand_search.best_estimator_
    #print(rand_search.best_params_)

    # creating and fitting the model created w/ optimized hyperparameters
    rf = RandomForestClassifier(max_depth=6, n_estimators=304, random_state=0)
    rf.fit(X_train, y_train)

    # testing the model 
    y_pred = rf.predict(X_test)

    # evaluating accuracy
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    print("Model Evaluation")
    print(f"Model accuracy: {accuracy}")
    print(f"Model precision: {precision}")
    print(f"Model recall: {recall}\n")

    #this fig is found in the ml_prediction/model_eval folder!
    #feature_importance = pd.Series(rf.feature_importances_, index=X_train.columns).sort_values(ascending=False)
    #feature_importance.plot.bar()
    #plt.title("Feature Importance")
    #plt.ylabel("Importance")
    #plt.xlabel("Features")
    #plt.show()

    # get actual result!!
    result = predict_from_smiles(args.smiles, rf)
    print(f"For {args.smiles}:\n")
    print(f"Molecular weight: {round(result[2][0], 3)}")
    print(f"Number of H-bond Acceptors: {result[2][1]}")
    print(f"Number of H-bond Donors: {result[2][2]}")
    print(f"Lipophilicty (LogP): {round(result[2][3],3)}")
    print("-------------------------------")
    print(f"Activity Prediction: {result[0]}")
    print(f"Probability: {result[1]}")

    

if __name__ == "__main__":
    main()


