import json

import model_ag
import model_ant


def _summarize(name: str, artifacts: dict):
    print(f"{name.upper()} MODEL")
    print("-" * (len(name) + 6))
    print("Internal holdout")
    print(json.dumps(artifacts["metrics"], indent=2, sort_keys=True))
    print("External scaffold validation")
    print(json.dumps(artifacts["external_validation"], indent=2, sort_keys=True))
    print("Benchmarks")
    print(json.dumps(artifacts["benchmarks"], indent=2, sort_keys=True))
    print()


def main():
    _summarize("agonist", model_ag.get_model_artifacts(force_retrain=True))
    _summarize("antagonist", model_ant.get_model_artifacts(force_retrain=True))


if __name__ == "__main__":
    main()
