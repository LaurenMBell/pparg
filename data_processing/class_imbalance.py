"""
I realized right about now that a set of jupyter notebooks 
would've worked just fine for this 

The split between active and inactive in the training data is 80% to 20%, and 
I want to deal with that. I had to lean pretty heavily on the tutorial for this part. 

Tutorials used:
- https://medium.com/data-science/class-imbalance-strategies-a-visual-guide-with-code-8bc8fae71e1a
"""

from sklearn.model_selection import train_test_split
import imblearn
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

#prepping the dataset
df = pd.read_csv("data/agonists.csv")
X_train, X_test, y_train, y_test = train_test_split(df.drop(["pChEMBL Value", "Standard Type"], axis = 1), 
                                                    df["pChEMBL Value"], 
                                                    test_size = 0.4, 
                                                    random_state=42)

#scale the dataset
scaler = MinMaxScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

#perform PCA to visualize the imbalanced dataset 
