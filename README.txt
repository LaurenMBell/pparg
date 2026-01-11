This is a side project of Lauren Bell's! It was inspired by a LinkedIn 
post I saw from Sumaiyyah Fatima* for the ML prediction pipeline. This is a QSAR-type pipeline, 
where users input a hypothesized drug and recieve prediction of that drug's effectiveness. 

To run the process, run pparg.py with your input molecule's SMILES code! 

Overview of the QSAR pipeline:
1. AQUISITION: Get PPARGs ChEMBL bioactivity dataset through web download:
    https://www.ebi.ac.uk/chembl/explore/activities/STATE_ID:TOB7olBkXsCNq8fsHg7p3g%3D%3D
2. PREPROCESSING: Clean it using Pandas to filter out missing data, duplicates,
    and label data according to bioactivity thresholds (to determine 'active'
    vs 'inactive' from IC50 value**). Then, remove salts from the interactions listed.
3. CONVERSION: Using RDKit, calculate molecular characteristics from the SMILES codes. Ex:
    molecular weight, lipophilicity, number of H-bond donors and acceptors.
4. EDA: Perform exploratory data analysis by comparing the distribution of interactions,
    and keep in mind Lipinski's rule of 5. Through correlation analysis, decide on relevant 
    features for the model. 
5. TRAINING: Using Scikit-learn, train a Random Forest model on the interaction data (one for 
    the agonists, and one for the antagonists)
6. EVALUATE/VALIDATE: Evaluate model performance using different metrics (F1, area under curve, confusion 
    matrix) to understand model precision
7. ML PREDICTION: The models should now be able to take a SMILES code as an input, and 
    output an activity prediction! 

------------------------------------------------------------------------------------------------------------
Notes:
* This is the LinkedIn post in question: 
https://www.linkedin.com/posts/sumaiyyahfatima_bioinformatics-machinelearning-drugdiscovery-activity-7351998459098398720-sskc?utm_source=share&utm_medium=member_desktop&rcm=ACoAAEj5d2sBvXjfVcYV5usp7owRrwpnnBEJ1L4

** I chose an IC50 (or, half-maximal inhibitory concentration) threshold of 10 uM, mostly 
informed by this forum from a 10 second google search lol: 
https://www.researchgate.net/post/Acceptable-IC50-drug-concentration-for-MTT-essay#:~:text=In%20most%20cases%2C%20the%20IC50,and%20excretion%20in%20our%20body.
I think it would be cool to play around with this threshold in a later iteration! 
And as a note, the value provided by ChEMBL is -log(molar IC50), so I just converted it from there to 5. 
So, pChEMBL >= 5 was active in my book, and pChEMBL < 5 was inactive
------------------------------------------------------------------------------------------------------------