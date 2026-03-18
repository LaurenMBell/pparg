import base64
import io
from dataclasses import dataclass
from typing import Any, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import requests
from flask import Flask, Response, render_template, request
from rdkit import Chem
from rdkit.Chem import Descriptors

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


def _is_probably_smiles(text: str) -> bool:
    if not text:
        return False
    m = Chem.MolFromSmiles(text)
    return m is not None


def pubchem_name_to_smiles(name: str) -> str:
    # PUG REST: SMILES for a compound name.
    # Note: PubChem may return ConnectivitySMILES instead of CanonicalSMILES
    # depending on dataset/availability, so we accept multiple fields.
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{requests.utils.quote(name)}/property/CanonicalSMILES/JSON"
    )
    r = requests.get(url, timeout=PUBCHEM_TIMEOUT_S)
    r.raise_for_status()
    data = r.json()
    props = data["PropertyTable"]["Properties"]
    if not props:
        raise ValueError("No PubChem results for that compound name.")
    row = props[0]
    if "CanonicalSMILES" in row:
        return row["CanonicalSMILES"]
    if "ConnectivitySMILES" in row:
        return row["ConnectivitySMILES"]
    if "IsomericSMILES" in row:
        return row["IsomericSMILES"]
    raise ValueError("PubChem response did not contain a SMILES field.")


def _draw_feature_importance(ag, ant):
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    ag_fi = ag["feature_importance"]
    ant_fi = ant["feature_importance"]
    features = list(ag_fi.index)

    x = np.arange(len(features))
    w = 0.36
    ax.bar(
        x - w / 2,
        [ag_fi[f] for f in features],
        width=w,
        color="orange",
        label="Agonist model",
    )
    ax.bar(
        x + w / 2,
        [ant_fi[f] for f in features],
        width=w,
        color="blue",
        label="Antagonist model",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(features, rotation=0)
    ax.set_ylabel("Importance")
    ax.set_title("Random Forest Feature Importance (training-derived)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper right")
    return _fig_to_data_uri(fig)


def _draw_model_metrics(ag, ant):
    fig, ax = plt.subplots(figsize=(7.6, 3.0))
    labels = ["accuracy", "precision", "recall"]
    ag_vals = [ag["metrics"][k] for k in labels]
    ant_vals = [ant["metrics"][k] for k in labels]
    x = np.arange(len(labels))
    w = 0.36
    ax.bar(x - w / 2, ag_vals, width=w, color="orange", label="Agonist model")
    ax.bar(x + w / 2, ant_vals, width=w, color="blue", label="Antagonist model")
    ax.set_xticks(x)
    ax.set_xticklabels([s.title() for s in labels])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Score")
    ax.set_title("Holdout Metrics (80/20 split, stratified, fixed seed)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper right")
    return _fig_to_data_uri(fig)


def _draw_prediction_hit_scores(ag_pred, ant_pred):
    fig, ax = plt.subplots(figsize=(7.6, 2.7))
    labels = ["Agonist", "Antagonist"]
    vals = [ag_pred["hit_score"], ant_pred["hit_score"]]
    colors = [c]
    ax.bar(labels, vals, color=colors)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Hit score (max class probability)")
    ax.set_title("Prediction Hit Scores")
    ax.grid(axis="y", alpha=0.25)
    for i, v in enumerate(vals):
        ax.text(i, min(0.98, v + 0.03), f"{v:.3f}", ha="center", va="bottom", fontsize=10)
    return _fig_to_data_uri(fig)


def _compute_chemical_traits(smiles: str) -> dict[str, Any]:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("Invalid SMILES: RDKit could not parse structure.")
    return {
        "heavy_atoms": int(m.GetNumHeavyAtoms()),
        "rings": int(Descriptors.RingCount(m)),
        "tpsa": float(Descriptors.TPSA(m)),
        "rot_bonds": int(Descriptors.NumRotatableBonds(m)),
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
    # One-sided: how often does the model assign a hit score >= this?
    ge = sum(1 for s in null_scores if s >= hit_score)
    return (ge + 1.0) / (len(null_scores) + 1.0)  # add-one smoothing

def _export_txt(
    compound_name: str,
    smiles: str,
    ag_pred: dict[str, Any],
    ant_pred: dict[str, Any],
    best: dict[str, Any],
    chem_traits: dict[str, Any],
    lipinski: dict[str, Any],
    stats: dict[str, Any],
) -> str:
    name_line = compound_name if compound_name else "—"
    lines = []
    lines.append("PPARG Activity Predictor — Export")
    lines.append("--------------------------------")
    lines.append(f"Compound name: {name_line}")
    lines.append(f"SMILES: {smiles}")
    lines.append("")
    lines.append("Predictions")
    lines.append(f"- Agonist: {ag_pred['activity']} (hit={ag_pred['hit_score']:.3f}, p={stats['ag_p_value']:.4f})")
    lines.append(
        f"- Antagonist: {ant_pred['activity']} (hit={ant_pred['hit_score']:.3f}, p={stats['ant_p_value']:.4f})"
    )
    lines.append(
        f"- Best hit: {best['model']} {best['activity']} (hit={best['hit_score']:.3f})"
    )
    lines.append("")
    lines.append("Descriptors (MW/HBA/HBD/LogP)")
    d = ag_pred["descriptors"]
    lines.append(f"- MW: {d['molecular_weight']:.2f}")
    lines.append(f"- HBA: {d['HBA']}")
    lines.append(f"- HBD: {d['HBD']}")
    lines.append(f"- LogP: {d['LogP']:.2f}")
    lines.append("")
    lines.append("Chemical traits")
    lines.append(f"- Heavy atoms: {chem_traits['heavy_atoms']}")
    lines.append(f"- Rings: {chem_traits['rings']}")
    lines.append(f"- tPSA: {chem_traits['tpsa']:.2f}")
    lines.append(f"- Rotatable bonds: {chem_traits['rot_bonds']}")
    lines.append("")
    lines.append("Lipinski checks")
    lines.append(f"- MW ≤ 500: {bool(lipinski['mw_le_500'])}")
    lines.append(f"- HBA ≤ 10: {bool(lipinski['hba_le_10'])}")
    lines.append(f"- HBD ≤ 5: {bool(lipinski['hbd_le_5'])}")
    lines.append(f"- LogP ≤ 5: {bool(lipinski['logp_le_5'])}")
    lines.append("")
    lines.append("Statistical framework (empirical)")
    lines.append(
        f"- Agonist null: n={stats['ag_null_n']} hit-scores from training distribution; "
        f"majority-baseline-acc={stats['ag_majority_acc']:.3f}"
    )
    lines.append(
        f"- Antagonist null: n={stats['ant_null_n']} hit-scores from training distribution; "
        f"majority-baseline-acc={stats['ant_majority_acc']:.3f}"
    )
    lines.append(
        "Note: p-values here are one-sided empirical scores for how extreme the model's hit score is "
        "relative to its own training distribution (not a biological assay significance test)."
    )
    lines.append("")
    return "\n".join(lines) + "\n"


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

        # choose "best hit" by probability; break ties by preferring active
        best = max(
            [ag_pred, ant_pred],
            key=lambda r: (r["hit_score"], 1 if r["activity"] == "active" else 0),
        )

        # chemical traits and lipinski flags from (shared) descriptors
        chem_traits = _compute_chemical_traits(smiles)
        lip = _lipinski_flags(ag_pred["descriptors"])

        # model artifacts for plots (cached)
        ag_art = model_ag.get_model_artifacts()
        ant_art = model_ant.get_model_artifacts()

        stats = {
            "ag_p_value": _empirical_p_value(ag_pred["hit_score"], ag_art["null"]["hit_scores"]),
            "ant_p_value": _empirical_p_value(ant_pred["hit_score"], ant_art["null"]["hit_scores"]),
            "ag_null_n": int(ag_art["null"]["n_train"]),
            "ant_null_n": int(ant_art["null"]["n_train"]),
            "ag_majority_acc": float(ag_art["null"]["majority_baseline_accuracy"]),
            "ant_majority_acc": float(ant_art["null"]["majority_baseline_accuracy"]),
        }

        figs = {
            "metrics": _draw_model_metrics(ag_art, ant_art),
            "feature_importance": _draw_feature_importance(ag_art, ant_art),
            "hit_scores": _draw_prediction_hit_scores(ag_pred, ant_pred),
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
        )
    except Exception as e:
        error = str(e)
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
        key=lambda r: (r["hit_score"], 1 if r["activity"] == "active" else 0),
    )
    chem_traits = _compute_chemical_traits(smiles)
    lip = _lipinski_flags(ag_pred["descriptors"])

    ag_art = model_ag.get_model_artifacts()
    ant_art = model_ant.get_model_artifacts()
    stats = {
        "ag_p_value": _empirical_p_value(ag_pred["hit_score"], ag_art["null"]["hit_scores"]),
        "ant_p_value": _empirical_p_value(ant_pred["hit_score"], ant_art["null"]["hit_scores"]),
        "ag_null_n": int(ag_art["null"]["n_train"]),
        "ant_null_n": int(ant_art["null"]["n_train"]),
        "ag_majority_acc": float(ag_art["null"]["majority_baseline_accuracy"]),
        "ant_majority_acc": float(ant_art["null"]["majority_baseline_accuracy"]),
    }

    txt = _export_txt(compound_name, smiles, ag_pred, ant_pred, best, chem_traits, lip, stats)
    filename = "pparg_results.txt"
    return Response(
        txt,
        mimetype="text/plain",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


if __name__ == "__main__":
    import os

    port = int(os.environ.get("PORT", "5051"))
    debug = os.environ.get("FLASK_DEBUG", "0") == "1"
    app.run(debug=debug, host="0.0.0.0", port=port)

