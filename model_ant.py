import argparse

from ml_core import ActivityModel


_MODEL = ActivityModel(model_name="antagonist", csv_path="data/antagonists.csv")


def get_model_artifacts(force_retrain: bool = False):
    return _MODEL.get_artifacts(force_retrain=force_retrain)


def predict(smiles: str):
    return _MODEL.predict(smiles)


def run_model_ant(smiles: str):
    result = predict(smiles)
    print(f"ANTAGONISTS: For {smiles}:")
    print("-------------------------------")
    print(f"Activity Prediction: {result['activity']}")
    print(f"Hit score: {result['hit_score']:.3f}")
    print(f"Calibrated P(active): {result['prob_active']:.3f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", required=True)
    args = parser.parse_args()

    artifacts = get_model_artifacts(force_retrain=True)
    metrics = artifacts["metrics"]
    external = artifacts["external_validation"]
    print("Antagonist Model Evaluation")
    print(f"Holdout accuracy: {metrics['accuracy']:.3f}")
    print(f"Holdout balanced accuracy: {metrics['balanced_accuracy']:.3f}")
    print(f"Holdout ROC-AUC: {metrics['roc_auc']:.3f}")
    print(f"External balanced accuracy: {external['balanced_accuracy']:.3f}")
    print(f"External ROC-AUC: {external['roc_auc']:.3f}\n")

    result = predict(args.smiles)
    print(f"For {args.smiles}:\n")
    for key, value in result["descriptors"].items():
        print(f"{key}: {value}")
    print("-------------------------------")
    print(f"Activity Prediction: {result['activity']}")
    print(f"Hit score: {result['hit_score']:.3f}")
    print(f"Calibrated P(active): {result['prob_active']:.3f}")


if __name__ == "__main__":
    main()
