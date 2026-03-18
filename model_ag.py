# Random Forest model to predict the EC50 score for agonists!

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, confusion_matrix
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from rdkit import Chem
from scipy.stats import randint
from rdkit.Chem import Descriptors, Lipinski, Crippen

import argparse

# -----------------------------
# Public, web-friendly API
# -----------------------------
_CACHE = {
    "rf": None,
    "metrics": None,
    "cm": None,
    "feature_importance": None,
    "null": None,
}

def _load_training_data():
    df = pd.read_csv("data/agonists.csv")
    df["activity"] = df["activity"].map({"inactive": 0, "active": 1})
    X = df[["molecular_weight", "HBA", "HBD", "LogP"]]
    y = df["activity"]
    return X, y

def get_model_artifacts(force_retrain: bool = False):
    """
    Returns cached artifacts (model + evaluation) for the agonist classifier.
    """
    if (not force_retrain) and _CACHE["rf"] is not None:
        return _CACHE

    X, y = _load_training_data()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    rf = RandomForestClassifier(max_depth=6, n_estimators=304, random_state=0)
    rf.fit(X_train, y_train)

    y_pred = rf.predict(X_test)
    metrics = {
        "accuracy": float(accuracy_score(y_test, y_pred)),
        "precision": float(precision_score(y_test, y_pred, zero_division=0)),
        "recall": float(recall_score(y_test, y_pred, zero_division=0)),
    }
    cm = confusion_matrix(y_test, y_pred, labels=[0, 1])
    feature_importance = pd.Series(
        rf.feature_importances_,
        index=X_train.columns,
    ).sort_values(ascending=False)

    # Empirical null for "hit score" (max predicted probability) and a simple null classifier.
    # This is NOT a mechanistic PPARG significance test; it quantifies how extreme a hit score is
    # relative to scores the model assigns on its own training distribution.
    p_train = rf.predict_proba(X)
    null_hit_scores = p_train.max(axis=1).astype(float).tolist()
    majority_class = int(y.value_counts().idxmax())
    null_classifier_accuracy = float((y_test == majority_class).mean())

    _CACHE.update(
        {
            "rf": rf,
            "metrics": metrics,
            "cm": cm,
            "feature_importance": feature_importance,
            "null": {
                "hit_scores": null_hit_scores,
                "majority_class": majority_class,
                "majority_baseline_accuracy": null_classifier_accuracy,
                "n_train": int(len(X)),
            },
        }
    )
    return _CACHE

def predict(smiles: str):
    """
    Predict agonist activity for a SMILES string.
    Returns a structured dict for UI use.
    """
    artifacts = get_model_artifacts()
    rf = artifacts["rf"]

    activity_label, hit_score, descriptors = predict_from_smiles(smiles, rf)
    return {
        "model": "agonist",
        "activity": activity_label,
        "hit_score": float(hit_score),
        "descriptors": {
            "molecular_weight": float(descriptors[0]),
            "HBA": int(descriptors[1]),
            "HBD": int(descriptors[2]),
            "LogP": float(descriptors[3]),
        },
        "metrics": artifacts["metrics"],
        "confusion_matrix": artifacts["cm"],
        "feature_importance": artifacts["feature_importance"],
    }

# same SMILES to descriptors conversion
def calc_descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("Invalid SMILES: RDKit could not parse structure.")
    
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

def run_model_ag(smiles):
     # reading in and mapping 
    artifacts = get_model_artifacts()
    result = predict_from_smiles(smiles, artifacts["rf"])
    print(f"AGONISTS: For {smiles}:")
    print("-------------------------------")
    print(f"Activity Prediction: {result[0]}")
    print(f"Probability: {result[1]}\n\n")
#===========================================================================

def main():
    # set up input args for the bash script
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles")

    args = parser.parse_args()

    artifacts = get_model_artifacts(force_retrain=True)
    metrics = artifacts["metrics"]
    print("Model Evaluation")
    print(f"Model accuracy: {metrics['accuracy']}")
    print(f"Model precision: {metrics['precision']}")
    print(f"Model recall: {metrics['recall']}\n")

    artifacts["feature_importance"].plot.bar()
    plt.title("Feature Importance")
    plt.ylabel("Importance")
    plt.xlabel("Features")
    plt.show()

    # get actual result!!
    result = predict_from_smiles(args.smiles, artifacts["rf"])
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


