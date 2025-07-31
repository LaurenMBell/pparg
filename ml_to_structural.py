# Get the hits from the ML models, decide molecular docking position, and pass to the structural side!!

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ag_pred")
    parser.add_argument("--ant_pred")

    args = parser.parse_args()

    if args.ag_pred == "active" and args.ant_pred == "inactive":
        type = "activator"
    elif args.ag_pred == "inactive" and args.ant_pred == "active":
        type = "inhibitor"
    elif args.ag_pred == "inactive" and args.ant_pred == "inactive":
        type = "ineffective"
    else:
        type = "mixed"
    
    

    

    

    
