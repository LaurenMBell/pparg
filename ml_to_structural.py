# Get the hits from the ML models, decide molecular docking position, and pass to the structural side!!

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ag_pred")
    parser.add_argument("--ant_pred")

    args = parser.parse_args()

    # 0 for activator/agonist, 1 for inhibitor/antagonist, 2 for ineffective
    # compound, 3 for mixed activity
    type = 0 

    if args.ag_pred == "active" and args.ant_pred == "inactive":
        type = 0
    elif args.ag_pred == "inactive" and args.ant_pred == "active":
        type = 1
    elif args.ag_pred == "inactive" and args.ant_pred == "inactive":
        type = 2
    else:
        type = 3

    

    
