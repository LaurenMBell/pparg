"""
tutorials used:
 - https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html
"""

import joblib 
import pandas as pd 

ag_model = joblib.load("agonist_rf_model.joblib")
ant_model = joblib.load("antagonist_rf_model.joblib")




