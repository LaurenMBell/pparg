"""
Tutorials/articles I used: 
 - https://towardsdatascience.com/fit-vs-predict-vs-fit-predict-in-python-scikit-learn-f15a34a8d39f/#:~:text=Note%20that%20every%20estimator%20might,you%20will%20receive%20a%20%5Bexceptions.
 - https://stackoverflow.com/questions/62646058/how-does-the-predict-method-work-on-scikit-learn
 - https://machinelearningmastery.com/make-predictions-scikit-learn/ 

for the time being, I think I'll build out a simpler CLI
and wrap a Flask app around it later
"""

import pandas as pd 
from data_processing.featSelection import Features
from sklearn.ensemble import RandomForestClassifier
import joblib
import numpy as np

def report(m, ag_hit, ant_hit, ag_prob, ant_prob):
    print("")
    if ag_hit==1 and ant_hit==1:
        if ag_prob > ant_prob: 
            print(f"{m} has predicted activity as a activator (ie. agonist).")
        else:
            print(f"{m} has predicted activity as a inhibitor (ie. antagonist).")

    if (ag_hit==0 and ant_hit==1):
        print(f"{m} has predicted activity as a inhibitor (ie. antagonist).")
    elif ag_hit==1 and ant_hit==0:
        print(f"{m} has predicted activity as a activator (ie. agonist).")
    elif ag_hit==0 and ant_hit==0:
        print(f"{m} is not predicted to have bioactivity with PPARG.")

    

def main():
    print("================= PPARG DRUG CANDIDATE PREDICTION =================")
    ag_rf = joblib.load("ml_prediction/agonist_rf_model.joblib")
    ant_rf = joblib.load("ml_prediction/antagonist_rf_model.joblib")

    print("Input SMILES to evaluate: ")
    smiles = input("> ")

    candidate_feats = pd.DataFrame([Features.calc_descriptors(smiles)]) #FEATURES RETURNS A DICT


    ag_pred = ag_rf.predict(candidate_feats)
    ant_pred = ant_rf.predict(candidate_feats)

    #probaility of class 1, active
    ag_probs = ag_rf.predict_proba(candidate_feats)[0][1]
    ant_probs = ant_rf.predict_proba(candidate_feats)[0][1]

    report(smiles, ag_pred, ant_pred, ag_probs, ant_probs)

    print("\nProbabilites (inactive, active):")
    print(f"Antagonist: {ant_rf.predict_proba(candidate_feats)[0]}")
    print(f"Agonist: {ag_rf.predict_proba(candidate_feats)[0]}")


if __name__ == "__main__":
    main()
