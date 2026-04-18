"""
Tutortials referenced: 
 - https://www.geeksforgeeks.org/machine-learning/random-forest-algorithm-in-machine-learning/

Until I build out a UI, a CLI will have to do :)
"""

import pandas as pd 
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

#undo what you did in class_imbalance.py 
X_train = pd.read_csv("../data_processing/data/agonists_training.csv")
