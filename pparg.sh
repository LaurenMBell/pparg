#!/bin/bash

# design actual pipeline logic here
# what do you need to pass to what part of the pipeline? 

<<PLAN
1)  bash will take SMILES as command line input, 
    and pass to both Ml models

    IN: SMILES code
    OUT: 2 activity predictions
    OUT: Evaluation metrics

2) bash will pass these predictions to ml_to_structual, 
    which will decide if its an inhibitor or an activator, 
    and spit out a molecular docking position

    IN: 2 activity predictions
    OUT: molecular docking position

3) bash will take the position and give it to the 
    structural side for docking 

    IN: molecular docking position
    OUT: structural prediction
    OUT: Evaluation metrics

4) bash will pass this, the activity predictions, and the evaluation
    metrics to final_report.py!
PLAN


