"""
Tutortials referenced: 
 - https://www.geeksforgeeks.org/machine-learning/random-forest-algorithm-in-machine-learning/

Until I build out a UI, a CLI will have to do :)
"""

import pandas as pd 
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

#undo what you did in class_imbalance.py 
X_train = pd.read_csv("../../data_processing/data/antagonists_training.csv")
y_train = X_train["activity"]
X_train = X_train.drop(columns=["activity"]) 

X_test = pd.read_csv("../../data_processing/data/antagonists_test.csv")
y_test = X_test["activity"]
X_test = X_test.drop(columns=["activity"])

ag_model = RandomForestClassifier(n_estimators=100,random_state=42)
ag_model.fit(X_train, y_train)

y_pred = ag_model.predict(X_test)

#model evaluation
accuracy = accuracy_score(y_test, y_pred)
classification = classification_report(y_test, y_pred)

print("Antagonist RF Model")
print(f"Accuracy: {accuracy:.2f}\n")
print(f"Classification: \n{classification}")
