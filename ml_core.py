from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
import os
import pickle
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.calibration import CalibratedClassifierCV
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_STATE = 42
FEATURE_COLUMNS = [
    "molecular_weight",
    "HBA",
    "HBD",
    "LogP",
    "heavy_atoms",
    "rings",
    "tpsa",
    "rot_bonds",
]
DISPLAY_FEATURE_NAMES = {
    "molecular_weight": "MW",
    "HBA": "HBA",
    "HBD": "HBD",
    "LogP": "LogP",
    "heavy_atoms": "Heavy atoms",
    "rings": "Rings",
    "tpsa": "tPSA",
    "rot_bonds": "Rotatable bonds",
}
RF_PARAMS = {
    "n_estimators": 200,
    "max_depth": 10,
    "min_samples_leaf": 2,
    "class_weight": "balanced_subsample",
    "random_state": RANDOM_STATE,
    "n_jobs": 1,
}


@dataclass(frozen=True)
class ValidationMetrics:
    accuracy: float
    balanced_accuracy: float
    precision: float
    recall: float
    f1: float
    roc_auc: float
    avg_precision: float
    brier_score: float


@dataclass(frozen=True)
class BenchmarkResult:
    model_name: str
    accuracy: float
    balanced_accuracy: float
    precision: float
    recall: float
    f1: float
    roc_auc: float
    avg_precision: float


class Featurizer(ABC):
    @abstractmethod
    def featurize_smiles(self, smiles: str) -> dict[str, float]:
        raise NotImplementedError

    def featurize_frame(self, smiles_series: pd.Series) -> pd.DataFrame:
        records = [self.featurize_smiles(smiles) for smiles in smiles_series]
        return pd.DataFrame(records, columns=FEATURE_COLUMNS)


class Predictor(ABC):
    @abstractmethod
    def predict(self, smiles: str) -> dict[str, Any]:
        raise NotImplementedError


class RDKitDescriptorFeaturizer(Featurizer):
    def featurize_smiles(self, smiles: str) -> dict[str, float]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES: RDKit could not parse structure.")

        return {
            "molecular_weight": float(Descriptors.MolWt(mol)),
            "HBA": float(Lipinski.NumHAcceptors(mol)),
            "HBD": float(Lipinski.NumHDonors(mol)),
            "LogP": float(Crippen.MolLogP(mol)),
            "heavy_atoms": float(mol.GetNumHeavyAtoms()),
            "rings": float(Descriptors.RingCount(mol)),
            "tpsa": float(Descriptors.TPSA(mol)),
            "rot_bonds": float(Descriptors.NumRotatableBonds(mol)),
        }


def metrics_to_dict(metrics: ValidationMetrics) -> dict[str, float]:
    return asdict(metrics)


def benchmark_to_dict(results: list[BenchmarkResult]) -> list[dict[str, float | str]]:
    return [asdict(result) for result in results]


def _compute_metrics(y_true: pd.Series, y_pred: np.ndarray, y_prob: np.ndarray) -> ValidationMetrics:
    return ValidationMetrics(
        accuracy=float(accuracy_score(y_true, y_pred)),
        balanced_accuracy=float(balanced_accuracy_score(y_true, y_pred)),
        precision=float(precision_score(y_true, y_pred, zero_division=0)),
        recall=float(recall_score(y_true, y_pred, zero_division=0)),
        f1=float(f1_score(y_true, y_pred, zero_division=0)),
        roc_auc=float(roc_auc_score(y_true, y_prob)),
        avg_precision=float(average_precision_score(y_true, y_prob)),
        brier_score=float(brier_score_loss(y_true, y_prob)),
    )


def _external_validation_mask(smiles: pd.Series, y: pd.Series, test_fraction: float = 0.2) -> pd.Series:
    holdout_mask = pd.Series(False, index=smiles.index)

    for label in sorted(y.unique()):
        idx = y[y == label].index
        label_smiles = smiles.loc[idx]
        scaffold_groups: dict[str, list[int]] = {}
        for row_idx, molecule_smiles in label_smiles.items():
            try:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(smiles=molecule_smiles)
            except Exception:
                scaffold = ""
            scaffold_groups.setdefault(scaffold or "__no_scaffold__", []).append(row_idx)

        target = max(1, int(round(len(idx) * test_fraction)))
        selected: list[int] = []
        for _, members in sorted(scaffold_groups.items(), key=lambda item: len(item[1]), reverse=True):
            if len(selected) >= target:
                break
            selected.extend(members)

        if len(selected) >= len(idx):
            fallback_n = max(1, min(len(idx) - 1, target))
            selected = list(idx[:fallback_n])

        holdout_mask.loc[selected] = True

    if holdout_mask.nunique() == 1:
        _, external_idx = train_test_split(
            smiles.index,
            test_size=test_fraction,
            random_state=RANDOM_STATE,
            stratify=y,
        )
        holdout_mask.loc[:] = False
        holdout_mask.loc[external_idx] = True

    return holdout_mask


def _build_benchmark_models() -> dict[str, Any]:
    return {
        "majority_baseline": DummyClassifier(strategy="most_frequent"),
        "logistic_balanced": Pipeline(
            steps=[
                ("scale", StandardScaler()),
                ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=RANDOM_STATE)),
            ]
        ),
        "random_forest_balanced": RandomForestClassifier(**RF_PARAMS),
    }


def _run_benchmarks(X: pd.DataFrame, y: pd.Series) -> list[BenchmarkResult]:
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)
    scoring = {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "precision": "precision",
        "recall": "recall",
        "f1": "f1",
        "roc_auc": "roc_auc",
        "avg_precision": "average_precision",
    }
    results: list[BenchmarkResult] = []
    for model_name, model in _build_benchmark_models().items():
        scores = cross_validate(model, X, y, cv=cv, scoring=scoring, n_jobs=1)
        results.append(
            BenchmarkResult(
                model_name=model_name,
                accuracy=float(scores["test_accuracy"].mean()),
                balanced_accuracy=float(scores["test_balanced_accuracy"].mean()),
                precision=float(scores["test_precision"].mean()),
                recall=float(scores["test_recall"].mean()),
                f1=float(scores["test_f1"].mean()),
                roc_auc=float(scores["test_roc_auc"].mean()),
                avg_precision=float(scores["test_avg_precision"].mean()),
            )
        )

    return sorted(results, key=lambda result: (result.balanced_accuracy, result.avg_precision), reverse=True)


def _majority_baseline_accuracy(y_true: pd.Series, training_labels: pd.Series) -> float:
    majority_class = int(training_labels.mode().iat[0])
    baseline = np.full(shape=len(y_true), fill_value=majority_class)
    return float(accuracy_score(y_true, baseline))


def load_activity_dataset(csv_path: str, featurizer: Featurizer | None = None) -> pd.DataFrame:
    featurizer = featurizer or RDKitDescriptorFeaturizer()
    df = pd.read_csv(csv_path)
    df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed:")]
    df = df.dropna(subset=["Smiles", "activity"]).copy()
    df["activity"] = df["activity"].map({"inactive": 0, "active": 1})
    df = df.dropna(subset=["activity"]).copy()
    feature_df = featurizer.featurize_frame(df["Smiles"])
    for column in FEATURE_COLUMNS:
        df[column] = feature_df[column].astype(float)
    df["activity"] = df["activity"].astype(int)
    return df


class ActivityModel(Predictor):
    def __init__(self, model_name: str, csv_path: str, force_retrain: bool = False):
        self.model_name = model_name
        self.csv_path = csv_path
        self.featurizer = RDKitDescriptorFeaturizer()
        self._artifacts: dict[str, Any] | None = None
        self.cache_path = os.path.join("model_cache", f"{model_name}_artifacts.pkl")
        if force_retrain:
            self.get_artifacts(force_retrain=True)

    def get_artifacts(self, force_retrain: bool = False) -> dict[str, Any]:
        if self._artifacts is not None and not force_retrain:
            return self._artifacts
        if not force_retrain and os.path.exists(self.cache_path):
            with open(self.cache_path, "rb") as handle:
                self._artifacts = pickle.load(handle)
            return self._artifacts

        df = load_activity_dataset(self.csv_path, self.featurizer)
        X = df[FEATURE_COLUMNS]
        y = df["activity"]

        external_mask = _external_validation_mask(df["Smiles"], y)
        dev_df = df.loc[~external_mask].copy()
        ext_df = df.loc[external_mask].copy()

        X_dev = dev_df[FEATURE_COLUMNS]
        y_dev = dev_df["activity"]
        X_ext = ext_df[FEATURE_COLUMNS]
        y_ext = ext_df["activity"]

        X_train, X_test, y_train, y_test = train_test_split(
            X_dev,
            y_dev,
            test_size=0.2,
            random_state=RANDOM_STATE,
            stratify=y_dev,
        )

        holdout_model = CalibratedClassifierCV(
            RandomForestClassifier(**RF_PARAMS),
            cv=3,
            method="sigmoid",
        )
        holdout_model.fit(X_train, y_train)
        y_test_prob = holdout_model.predict_proba(X_test)[:, 1]
        y_test_pred = holdout_model.predict(X_test)
        holdout_metrics = _compute_metrics(y_test, y_test_pred, y_test_prob)
        cm = confusion_matrix(y_test, y_test_pred, labels=[0, 1])

        external_model = CalibratedClassifierCV(
            RandomForestClassifier(**RF_PARAMS),
            cv=3,
            method="sigmoid",
        )
        external_model.fit(X_dev, y_dev)
        y_ext_prob = external_model.predict_proba(X_ext)[:, 1]
        y_ext_pred = external_model.predict(X_ext)
        external_metrics = _compute_metrics(y_ext, y_ext_pred, y_ext_prob)

        final_model = CalibratedClassifierCV(
            RandomForestClassifier(**RF_PARAMS),
            cv=3,
            method="sigmoid",
        )
        final_model.fit(X, y)

        feature_model = RandomForestClassifier(**RF_PARAMS)
        feature_model.fit(X, y)
        feature_importance = pd.Series(
            feature_model.feature_importances_,
            index=FEATURE_COLUMNS,
        ).sort_values(ascending=False)

        full_prob = final_model.predict_proba(X)[:, 1]
        null_hit_scores = np.maximum(full_prob, 1.0 - full_prob).astype(float).tolist()
        benchmarks = _run_benchmarks(X_dev, y_dev)

        self._artifacts = {
            "model": final_model,
            "metrics": metrics_to_dict(holdout_metrics),
            "external_validation": metrics_to_dict(external_metrics),
            "confusion_matrix": cm,
            "feature_importance": feature_importance,
            "feature_columns": FEATURE_COLUMNS,
            "benchmarks": benchmark_to_dict(benchmarks),
            "null": {
                "hit_scores": null_hit_scores,
                "majority_class": int(y.mode().iat[0]),
                "majority_baseline_accuracy": _majority_baseline_accuracy(y_test, y_train),
                "external_majority_baseline_accuracy": _majority_baseline_accuracy(y_ext, y_dev),
                "n_train": int(len(X)),
            },
            "splits": {
                "development_rows": int(len(dev_df)),
                "external_rows": int(len(ext_df)),
            },
        }
        os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
        with open(self.cache_path, "wb") as handle:
            pickle.dump(self._artifacts, handle)
        return self._artifacts

    def predict(self, smiles: str) -> dict[str, Any]:
        artifacts = self.get_artifacts()
        descriptors = self.featurizer.featurize_smiles(smiles)
        ddf = pd.DataFrame([[descriptors[column] for column in FEATURE_COLUMNS]], columns=FEATURE_COLUMNS)
        prob_active = float(artifacts["model"].predict_proba(ddf)[0, 1])
        prediction = int(prob_active >= 0.5)
        hit_score = prob_active if prediction == 1 else 1.0 - prob_active

        descriptor_payload: dict[str, Any] = {}
        for key, value in descriptors.items():
            if key in {"HBA", "HBD", "heavy_atoms", "rings", "rot_bonds"}:
                descriptor_payload[key] = int(round(value))
            else:
                descriptor_payload[key] = float(value)

        return {
            "model": self.model_name,
            "activity": "active" if prediction == 1 else "inactive",
            "hit_score": float(hit_score),
            "prob_active": prob_active,
            "descriptors": descriptor_payload,
            "metrics": artifacts["metrics"],
            "external_validation": artifacts["external_validation"],
            "confusion_matrix": artifacts["confusion_matrix"],
            "feature_importance": artifacts["feature_importance"],
            "benchmarks": artifacts["benchmarks"],
        }
