This is a practice project of Lauren Bell's! It was inspired by a LinkedIn 
post I saw from Sumaiyyah Fatima* for the ML prediction pipeline, but the 
structural pipeline was created by me (sort of, it seems like the standard
process for molecular docking simulation) and the two were integrated together 
into a final prediction by me! This is a QSAR-type pipeline, and I hope to create a 
web framework for users to input their own SMILES codes, and recieve a drug prediction
as a result as the next part of this project.

To run the process, run pparg.sh with your input molecule's SMILES code! 


Overview of both sides of the QSAR pipeline and their integration:

ML PIPELINE (adapted from Sumaiyyah Fatima's schema):

1. AQUISITION: Get PPARGs ChEMBL bioactivity dataset through web download
2. PREPROCESSING: Clean it using Pandas to filter out missing data, duplicates,
    and label data according to bioactivity thresholds (to determine 'active'
    vs 'inactive' from IC50 value**). Then, remove salts from the interactions listed.
3. CONVERSION: Using RDKit, calculate molecular characteristics from the SMILES codes. Ex:
    molecular weight, lipophilicity, number of H-bond donors and acceptors, topological
    surface area. 
4. EDA: Perform exploratory data analysis by comparing the distribution of interactions,
    and keep in mind Lipinski's rule of 5. Through correlation analysis, decide on relevant 
    features for the model. 
5. TRAINING: Using Scikit-learn, train a Random Forest model on the interaction data. 
6. EVALUATE/VALIDATE: Evaluate model performance using different metrics (F1, area under curve, confusion 
    matrix) to understand model precision
7. ML PREDICTION: The model should now be able to take a SMILES code as an input, and 
    output an activity prediction! 


STRUCTURAL PIPELINE: 

1. Using the AlphaFold PPARG structure, use the protein coordinates structure 
2. Perform molecular docking for each interaction between PPARG and the drug candidate's
    structure - receive a predicted binding pose and affinity score. 
3. Use this data (structural interaction predictions) in the final step of the pipeline


COMBINED PREDICTION + ANALYSIS = INTEGRATION!


------------------------------------------------------------------------------------------------------------
Notes:
* This is the LinkedIn post in question: 
https://www.linkedin.com/posts/sumaiyyahfatima_bioinformatics-machinelearning-drugdiscovery-activity-7351998459098398720-sskc?utm_source=share&utm_medium=member_desktop&rcm=ACoAAEj5d2sBvXjfVcYV5usp7owRrwpnnBEJ1L4

** I chose an IC50 (or, half-maximal inhibitory concentration) threshold of 10 uM, mostly 
informed by this forum from a 10 second google search lol: 
https://www.researchgate.net/post/Acceptable-IC50-drug-concentration-for-MTT-essay#:~:text=In%20most%20cases%2C%20the%20IC50,and%20excretion%20in%20our%20body.
I think it would be cool to play around with this threshold in a later iteration!
------------------------------------------------------------------------------------------------------------