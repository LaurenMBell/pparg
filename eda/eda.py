"""
4/16/26

Tutorials/books referenced: 
 - https://medium.com/@post.gourang/a-holistic-guide-to-exploratory-data-analysis-eda-for-machine-learning-and-deep-learning-bc4f18f0143b
 - https://www.linkedin.com/pulse/best-python-libraries-exploratory-data-analysis-eda-suraj-kumar-soni-pqruc/
 - http://ciml.info/dl/v0_99/ciml-v0_99-ch05.pdf
 
What is the problem I'm aiming to solve here? 
    I want to classify molecule-molecule interactions using random forests. Or, 
    I guess put more simply, I want to classify bioactivity (IC50) of untested
    compounds using the activity data of tested compounds. 
 """

import matplotlib.pyplot as plt
import pandas as pd 
import sweetviz as sv

def summary_stats(df):
    print(f"Number of (rows, cols): {df.shape}")
    
    report = sv.analyze([df, 'Train'], target_feat = "pChEMBL Value")
    report.show_html('Report.html')


def main():
    bioactivity = pd.read_csv("../preprocessing/data/bioactivity.tsv", delimiter="\t")
    bioactivity = bioactivity.dropna(subset="pChEMBL Value") #drop any NaN IC50/EC50 values

    summary_stats(bioactivity) #create EDA report with sweetviz
    return 

if __name__ == "__main__":
    main()