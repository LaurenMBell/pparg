"""
Tutortials/docs referenced: 
 - https://www.geeksforgeeks.org/machine-learning/random-forest-algorithm-in-machine-learning/
 - https://joblib.readthedocs.io/en/stable/
 - https://www.analyticsvidhya.com/blog/2023/02/how-to-save-and-load-machine-learning-models-in-python-using-joblib-library/


Until I build out a UI, a CLI will have to do :)

Using JobLib to save the model 
"""

import pandas as pd 
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import joblib

class AntagonistRF:
    def __init__(self):
        #load data when calss is initialized instead of retraining with each call
        #undo what you did in class_imbalance.py 

        self.X_train = pd.read_csv("../../data_processing/data/antagonists_training.csv")
        self.y_train = self.X_train["activity"]
        self.X_train = self.X_train.drop(columns=["activity"]) 

        self.X_test = pd.read_csv("../../data_processing/data/antagonists_test.csv")
        self.y_test = self.X_test["activity"]
        self.X_test = self.X_test.drop(columns=["activity"])

        self.ag_model = RandomForestClassifier(n_estimators=100,random_state=42)
        self.ag_model.fit(self.X_train, self.y_train)

    def train(self):
        self.ag_model.fit(self.X_train, self.y_train)
    
    def basic_eval(self):
        y_pred = self.ag_model.predict(self.X_test)

        #model evaluation
        accuracy = accuracy_score(self.y_test, y_pred)
        classification = classification_report(self.y_test, y_pred)

        print("Antagonist RF Model")
        print(f"Accuracy: {accuracy:.2f}\n")
        print(f"Classification: \n{classification}")

    def save_model(self):
        joblib.dump(self.ag_model, "antagonist_rf_model.joblib")



