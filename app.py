import base64
import io
import os
from typing import Any, Optional

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import requests
from flask import Flask, Response, render_template, request
from matplotlib import colors as mcolors
from rdkit import Chem

from ml_core import DISPLAY_FEATURE_NAMES
import model_ag
import model_ant


app = Flask(__name__)


PUBCHEM_TIMEOUT_S = 12


def _fig_to_data_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def _paired_oranges() -> tuple[str, str]:
    cmap = plt.get_cmap("Paired")
    light = mcolors.to_hex(cmap(6))
    dark = mcolors.to_hex(cmap(7))
    return light, dark


def _is_probably_smiles(text: str) -> bool:
    if not text:
        return False
    return Chem.MolFromSmiles(text) is not None


def pubchem_name_to_smiles(name: str) -> str:
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{requests.utils.quote(name)}/property/CanonicalSMILES/JSON"
    )
    response = requests.get(url, timeout=PUBCHEM_TIMEOUT_S)
    response.raise_for_status()
    data = response.json()
    props = data["PropertyTable"]["Properties"]
    if not props:
        raise ValueError("No PubChem results for that compound name.")
    row = props[0]
    for key in ("CanonicalSMILES", "ConnectivitySMILES", "IsomericSMILES"):
        if key in row:
            return row[key]
    raise ValueError("PubChem response did not contain a SMILES field.")


def _draw_feature_importance(ag, ant):
    orange_light, orange_dark = _paired_oranges()
    fig, ax = plt.subplots(figsize=(9.4, 3.6))
    ag_fi = ag["feature_importance"]
    ant_fi = ant["feature_importance"]
    features = [DISPLAY_FEATURE_NAMES.get(feature, feature) for feature in ag_fi.index]

    x = np.arange(len(features))
    w = 0.36
    ax.bar(x - w / 2, [ag_fi[f] for f in ag_fi.index], width=w, color=orange_dark, label="Agonist model")
    ax.bar(x + w / 2, [ant_fi[f] for f in ant_fi.index], width=w, color=orange_light, label="Antagonist model")

    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=0)
    ax.set_ylabel("Importance")
    ax.set_title("Random Forest Feature Importance (training-derived)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper right")
    return _fig_to_data_uri(fig)


def _draw_model_metrics(ag, ant):
    orange_light, orange_dark = _paired_oranges()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.2), constrained_layout=True)
    labels = ["balanced_accuracy", "roc_auc", "avg_precision"]
    pretty = ["Balanced Acc.", "ROC-AUC", "PR-AUC"]
    x = np.arange(len(labels))
    w = 0.36

    for ax, payload, title in (
        (axes[0], ag, "Agonist validation"),
        (axes[1], ant, "Antagonist validation"),
    ):
        internal = [payload["metrics"][key] for key in labels]
        external = [payload["external_validation"][key] for key in labels]
        ax.bar(x - w / 2, internal, width=w, color=orange_dark, label="Internal holdout")
        ax.bar(x + w / 2, external, width=w, color=orange_light, label="External scaffold")
        ax.set_xticks(x)
        ax.set_xticklabels(pretty)
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("Score")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25)
        ax.legend(frameon=False, loc="lower right")

    return _fig_to_data_uri(fig)


def _draw_prediction_hit_scores(ag_pred, ant_pred):
    orange_light, orange_dark = _paired_oranges()
    fig, ax = plt.subplots(figsize=(7.6, 2.7))
    labels = ["Agonist", "Antagonist"]
    vals = [ag_pred["hit_score"], ant_pred["hit_score"]]
    ax.bar(labels, vals, color=[orange_dark, orange_light])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Hit score")
    ax.set_title("Prediction Confidence")
    ax.grid(axis="y", alpha=0.25)
    for i, value in enumerate(vals):
        ax.text(i, min(0.98, value + 0.03), f"{value:.3f}", ha="center", va="bottom", fontsize=10)
    return _fig_to_data_uri(fig)


def _draw_benchmark_summary(ag, ant):
    orange_light, orange_dark = _paired_oranges()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.4), constrained_layout=True)

    for ax, payload, title in (
        (axes[0], ag, "Agonist benchmark"),
        (axes[1], ant, "Antagonist benchmark"),
    ):
        bench = payload["benchmarks"]
        labels = [row["model_name"].replace("_", "\n") for row in bench]
        vals = [row["balanced_accuracy"] for row in bench]
        ax.bar(labels, vals, color=[orange_dark, orange_light, "#b6b6b6"][: len(vals)])
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("Balanced accuracy")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25)

    return _fig_to_data_uri(fig)


def _descriptor_items(descriptors: dict[str, Any]) -> list[tuple[str, Any]]:
    rows = []
    for key, label in DISPLAY_FEATURE_NAMES.items():
        value = descriptors[key]
        rows.append((label, round(value, 2) if isinstance(value, float) else value))
    return rows


def _compute_chemical_traits(descriptors: dict[str, Any]) -> dict[str, Any]:
    return {
        "heavy_atoms": int(descriptors["heavy_atoms"]),
        "rings": int(descriptors["rings"]),
        "tpsa": float(descriptors["tpsa"]),
        "rot_bonds": int(descriptors["rot_bonds"]),
    }


def _lipinski_flags(descriptors: dict[str, Any]) -> dict[str, bool]:
    return {
        "mw_le_500": descriptors["molecular_weight"] <= 500,
        "hba_le_10": descriptors["HBA"] <= 10,
        "hbd_le_5": descriptors["HBD"] <= 5,
        "logp_le_5": descriptors["LogP"] <= 5,
    }


def _empirical_p_value(hit_score: float, null_scores: list[float]) -> float:
    if not null_scores:
        return float("nan")
    ge = sum(1 for score in null_scores if score >= hit_score)
    return (ge + 1.0) / (len(null_scores) + 1.0)


def _validation_snapshot(artifacts: dict[str, Any]) -> dict[str, Any]:
    benchmarks = artifacts["benchmarks"]
    return {
        "internal": artifacts["metrics"],
        "external": artifacts["external_validation"],
        "benchmarks": benchmarks,
        "best_benchmark": benchmarks[0],
        "splits": artifacts["splits"],
    }


def _export_txt(
    compound_name: str,
    smiles: str,
    ag_pred: dict[str, Any],
    ant_pred: dict[str, Any],
    best: dict[str, Any],
    chem_traits: dict[str, Any],
    lipinski: dict[str, Any],
    stats: dict[str, Any],
    validation: dict[str, Any],
) -> str:
    name_line = compound_name if compound_name else "—"
    lines = [
        "PPARG Activity Predictor — Export",
        "--------------------------------",
        f"Compound name: {name_line}",
        f"SMILES: {smiles}",
        "",
        "Predictions",
        f"- Agonist: {ag_pred['activity']} (hit={ag_pred['hit_score']:.3f}, p={stats['ag_p_value']:.4f})",
        f"- Antagonist: {ant_pred['activity']} (hit={ant_pred['hit_score']:.3f}, p={stats['ant_p_value']:.4f})",
        f"- Best hit: {best['model']} {best['activity']} (hit={best['hit_score']:.3f})",
        "",
        "Descriptors",
    ]
    for label, value in _descriptor_items(ag_pred["descriptors"]):
        if isinstance(value, float):
            lines.append(f"- {label}: {value:.2f}")
        else:
            lines.append(f"- {label}: {value}")

    lines.extend(
        [
            "",
            "Chemical traits",
            f"- Heavy atoms: {chem_traits['heavy_atoms']}",
            f"- Rings: {chem_traits['rings']}",
            f"- tPSA: {chem_traits['tpsa']:.2f}",
            f"- Rotatable bonds: {chem_traits['rot_bonds']}",
            "",
            "Validation",
            (
                f"- Internal holdout: balanced_acc={validation['internal']['balanced_accuracy']:.3f}, "
                f"roc_auc={validation['internal']['roc_auc']:.3f}, "
                f"pr_auc={validation['internal']['avg_precision']:.3f}"
            ),
            (
                f"- External scaffold split: balanced_acc={validation['external']['balanced_accuracy']:.3f}, "
                f"roc_auc={validation['external']['roc_auc']:.3f}, "
                f"pr_auc={validation['external']['avg_precision']:.3f}"
            ),
            (
                f"- Best benchmark model: {validation['best_benchmark']['model_name']} "
                f"(balanced_acc={validation['best_benchmark']['balanced_accuracy']:.3f})"
            ),
            "",
            "Lipinski checks",
            f"- MW ≤ 500: {bool(lipinski['mw_le_500'])}",
            f"- HBA ≤ 10: {bool(lipinski['hba_le_10'])}",
            f"- HBD ≤ 5: {bool(lipinski['hbd_le_5'])}",
            f"- LogP ≤ 5: {bool(lipinski['logp_le_5'])}",
            "",
            "Statistical framework (empirical)",
            (
                f"- Agonist null: n={stats['ag_null_n']} hit-scores from training distribution; "
                f"internal-majority-acc={stats['ag_majority_acc']:.3f}; "
                f"external-majority-acc={stats['ag_external_majority_acc']:.3f}"
            ),
            (
                f"- Antagonist null: n={stats['ant_null_n']} hit-scores from training distribution; "
                f"internal-majority-acc={stats['ant_majority_acc']:.3f}; "
                f"external-majority-acc={stats['ant_external_majority_acc']:.3f}"
            ),
            (
                "Note: p-values here are one-sided empirical scores for how extreme the model's hit score is "
                "relative to its own training distribution (not a biological assay significance test)."
            ),
            "",
        ]
    )

    return "\n".join(lines)


@app.get("/")
def index():
    return render_template("index.html")


@app.post("/predict")
def predict():
    compound_name = (request.form.get("compound_name") or "").strip()
    smiles_in = (request.form.get("smiles") or "").strip()

    error: Optional[str] = None
    smiles: Optional[str] = None

    try:
        if smiles_in:
            if not _is_probably_smiles(smiles_in):
                raise ValueError("Provided SMILES could not be parsed.")
            smiles = smiles_in
        else:
            if not compound_name:
                raise ValueError("Provide either a compound name or a SMILES string.")
            smiles = pubchem_name_to_smiles(compound_name)

        ag_pred = model_ag.predict(smiles)
        ant_pred = model_ant.predict(smiles)
        best = max(
            [ag_pred, ant_pred],
            key=lambda result: (result["hit_score"], 1 if result["activity"] == "active" else 0),
        )

        chem_traits = _compute_chemical_traits(ag_pred["descriptors"])
        lip = _lipinski_flags(ag_pred["descriptors"])

        ag_art = model_ag.get_model_artifacts()
        ant_art = model_ant.get_model_artifacts()
        ag_validation = _validation_snapshot(ag_art)
        ant_validation = _validation_snapshot(ant_art)

        stats = {
            "ag_p_value": _empirical_p_value(ag_pred["hit_score"], ag_art["null"]["hit_scores"]),
            "ant_p_value": _empirical_p_value(ant_pred["hit_score"], ant_art["null"]["hit_scores"]),
            "ag_null_n": int(ag_art["null"]["n_train"]),
            "ant_null_n": int(ant_art["null"]["n_train"]),
            "ag_majority_acc": float(ag_art["null"]["majority_baseline_accuracy"]),
            "ant_majority_acc": float(ant_art["null"]["majority_baseline_accuracy"]),
            "ag_external_majority_acc": float(ag_art["null"]["external_majority_baseline_accuracy"]),
            "ant_external_majority_acc": float(ant_art["null"]["external_majority_baseline_accuracy"]),
        }

        figs = {
            "metrics": _draw_model_metrics(ag_art, ant_art),
            "feature_importance": _draw_feature_importance(ag_art, ant_art),
            "hit_scores": _draw_prediction_hit_scores(ag_pred, ant_pred),
            "benchmarks": _draw_benchmark_summary(ag_art, ant_art),
        }

        return render_template(
            "results.html",
            compound_name=compound_name,
            smiles=smiles,
            ag=ag_pred,
            ant=ant_pred,
            best=best,
            chem_traits=chem_traits,
            lipinski=lip,
            figs=figs,
            stats=stats,
            ag_validation=ag_validation,
            ant_validation=ant_validation,
            descriptor_items=_descriptor_items(ag_pred["descriptors"]),
        )
    except Exception as exc:
        error = str(exc)
        return render_template(
            "index.html",
            error=error,
            compound_name=compound_name,
            smiles=smiles_in,
        )


@app.post("/export")
def export():
    compound_name = (request.form.get("compound_name") or "").strip()
    smiles_in = (request.form.get("smiles") or "").strip()
    if not smiles_in:
        return Response("Missing SMILES.\n", status=400, mimetype="text/plain")

    smiles = smiles_in
    ag_pred = model_ag.predict(smiles)
    ant_pred = model_ant.predict(smiles)
    best = max(
        [ag_pred, ant_pred],
        key=lambda result: (result["hit_score"], 1 if result["activity"] == "active" else 0),
    )
    chem_traits = _compute_chemical_traits(ag_pred["descriptors"])
    lip = _lipinski_flags(ag_pred["descriptors"])

    ag_art = model_ag.get_model_artifacts()
    ant_art = model_ant.get_model_artifacts()
    validation = _validation_snapshot(ag_art)
    stats = {
        "ag_p_value": _empirical_p_value(ag_pred["hit_score"], ag_art["null"]["hit_scores"]),
        "ant_p_value": _empirical_p_value(ant_pred["hit_score"], ant_art["null"]["hit_scores"]),
        "ag_null_n": int(ag_art["null"]["n_train"]),
        "ant_null_n": int(ant_art["null"]["n_train"]),
        "ag_majority_acc": float(ag_art["null"]["majority_baseline_accuracy"]),
        "ant_majority_acc": float(ant_art["null"]["majority_baseline_accuracy"]),
        "ag_external_majority_acc": float(ag_art["null"]["external_majority_baseline_accuracy"]),
        "ant_external_majority_acc": float(ant_art["null"]["external_majority_baseline_accuracy"]),
    }

    txt = _export_txt(compound_name, smiles, ag_pred, ant_pred, best, chem_traits, lip, stats, validation)
    return Response(
        txt,
        mimetype="text/plain",
        headers={"Content-Disposition": 'attachment; filename="pparg_results.txt"'},
    )


if __name__ == "__main__":
    import os

    port = int(os.environ.get("PORT", "5051"))
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(debug=debug, host="0.0.0.0", port=port)
