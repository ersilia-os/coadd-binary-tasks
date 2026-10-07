"""
Validate CoADD models against EU-OPENSCREEN experimental data.

Step 11 predicted the EU-OPENSCREEN library with every inhib_50 / mic_25 model.
EU-OPENSCREEN also ships *measured* activity per pathogen (active/inactive), so
this step closes the loop: for each model, compare its predictions to the
experimental labels and quantify agreement with AUROC.

For each step-11 prediction file:
  1. Attach the experimental label by merging predictions with the per-pathogen
     EU-OPENSCREEN file on SMILES (also brings in the InChIKey).
  2. Exclude any compound present in *that model's own* training set, matched by
     InChIKey (EU-OPENSCREEN SMILES are a different string format than the
     training std_smiles, so raw-SMILES matching would miss them).
  3. Drop ambiguous labels (bin == -1).
  4. Compute AUROC of y_hat (predict_proba) vs the binary experimental label.

Output: output/12_validate_euopenscreen/
  12_auroc_results.csv   one row per model
  12_roc_curves.png      one ROC panel per pathogen (a curve per model)
  12_auroc_barplot.png   AUROC per model, grouped by pathogen

Run with:
  conda run -n h3d python scripts/12_validate_euopenscreen.py
"""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd
import stylia as st
from matplotlib.patches import Patch
from sklearn.metrics import roc_auc_score, roc_curve

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

# Format: print | Style: ersilia
st.set_format("print")
st.set_style("ersilia")
nc = st.NamedColors()

# One colour per measure, consistent with script 11.
MEASURE_COLORS = {"inhib": nc.blue, "mic": nc.orange}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
euos_dir   = os.path.join(root, "..", "data", "raw", "eu-openscreen")
coadd_dir  = os.path.join(root, "..", "data", "processed", "coadd")
preds_dir  = os.path.join(root, "..", "output", "11_predict_euopenscreen")
output_dir = os.path.join(root, "..", "output", "12_validate_euopenscreen")
os.makedirs(output_dir, exist_ok=True)

ASSAY_DIR = {"inhib": "03_binarised_inhibition", "mic": "05_binarised_mic"}


def parse_model_key(name):
    """Split 'patho_strain_assay_cutoff' into components."""
    parts = name.rsplit("_", 1)
    cutoff = parts[1]
    rest   = parts[0].rsplit("_", 1)
    assay  = rest[1]
    rest2  = rest[0].split("_", 1)
    patho  = rest2[0]
    strain = rest2[1]
    return patho, strain, assay, cutoff


def training_inchikeys(patho, strain, assay):
    """InChIKeys of the compounds this model was trained on."""
    fpath = os.path.join(coadd_dir, ASSAY_DIR[assay], f"{patho}_{strain}.csv")
    if not os.path.exists(fpath):
        print(f"  WARNING: training file missing ({os.path.basename(fpath)})")
        return set()
    return set(pd.read_csv(fpath, usecols=["inchikey"])["inchikey"])


# ---------------------------------------------------------------------------
# Per-model validation
# ---------------------------------------------------------------------------
pred_files = sorted(glob.glob(os.path.join(preds_dir, "11_*_predictions.json")))
print(f"Prediction files: {len(pred_files)}")

# Cache experimental tables (one per pathogen) to avoid re-reading.
exp_cache = {}


def load_experimental(patho):
    if patho not in exp_cache:
        exp_cache[patho] = pd.read_csv(
            os.path.join(euos_dir, f"02_{patho}.csv")
        )[["smiles", "inchikey", "bin"]]
    return exp_cache[patho]


results = []        # rows for the summary CSV
roc_curves = {}     # model_key -> (fpr, tpr) for plotting

for pred_path in pred_files:
    model_key = os.path.basename(pred_path)[len("11_"):-len("_predictions.json")]
    patho, strain, assay, cutoff = parse_model_key(model_key)

    with open(pred_path) as fh:
        pred = json.load(fh)
    pred_df = pd.DataFrame({"smiles": pred["smiles"], "y_hat": pred["y_hat"]})

    # Attach experimental label + inchikey via SMILES (same source as step 11).
    exp = load_experimental(patho)
    df = pred_df.merge(exp, on="smiles", how="inner")

    # Exclude this model's own training compounds (by InChIKey).
    train_iks = training_inchikeys(patho, strain, assay)
    in_train = df["inchikey"].isin(train_iks)
    n_excluded = int(in_train.sum())
    df = df[~in_train]

    # Keep only definitive labels.
    df = df[df["bin"].isin([0, 1])]

    # Drop compounds the model could not score. lazyqsar returns NaN for SMILES RDKit
    # cannot parse rather than scoring them anyway; in this library that is 19 carborane
    # cages out of 106,317 (0.018%), the same ones for every model. roc_auc_score raises
    # on NaN, so they are removed here instead of poisoning the whole AUROC.
    n_unscored = int(df["y_hat"].isna().sum())
    if n_unscored:
        df = df[df["y_hat"].notna()]

    n_eval = len(df)
    n_active = int((df["bin"] == 1).sum())
    if df["bin"].nunique() < 2:
        print(f"[{model_key}] only one class after filtering — skipping.")
        continue

    auroc = roc_auc_score(df["bin"], df["y_hat"])
    fpr, tpr, _ = roc_curve(df["bin"], df["y_hat"])
    roc_curves[model_key] = (fpr, tpr)

    results.append({
        "n_unscored":       n_unscored,
        "model_key":        model_key,
        "patho":            patho,
        "strain":           strain,
        "assay":            assay,
        "cutoff":           cutoff,
        "n_eval":           n_eval,
        "n_active":         n_active,
        "n_excluded_train": n_excluded,
        "auroc":            round(auroc, 4),
    })
    print(f"[{model_key}] n_eval={n_eval:,} n_active={n_active} "
          f"excluded={n_excluded} AUROC={auroc:.4f}")

results_df = pd.DataFrame(results)
csv_path = os.path.join(output_dir, "12_auroc_results.csv")
results_df.to_csv(csv_path, index=False)
print(f"\nSaved: {os.path.basename(csv_path)}  ({len(results_df)} models)")

# ---------------------------------------------------------------------------
# Plot A: ROC curves, one panel per pathogen
# ---------------------------------------------------------------------------
pathos = sorted(results_df["patho"].unique())


def plot_roc_models(ax, rows):
    """One ROC curve per model in `rows` (a per-pathogen slice of results_df)."""
    pal = st.CategoricalPalette("ersilia")
    colors = pal.get(len(rows))
    ax.plot([0, 1], [0, 1], linestyle="--", color=nc.gray)
    for color, (_, r) in zip(colors, rows.iterrows()):
        fpr, tpr = roc_curves[r["model_key"]]
        ax.plot(fpr, tpr, color=color,
                label=f"{r['strain']}_{r['assay']} AUC={r['auroc']:.3f}")
    ax.legend(loc="lower right")


_, axs = st.create_figure(2, 3)
for patho in pathos:
    ax = axs.next()
    rows = results_df[results_df["patho"] == patho]
    plot_roc_models(ax, rows)
    st.label(ax, xlabel="FPR", ylabel="TPR", title=patho)
for _ in range(6 - len(pathos)):
    axs.next().axis("off")

roc_path = os.path.join(output_dir, "12_roc_curves.png")
st.save_figure(roc_path)
print(f"Saved: {os.path.basename(roc_path)}")

# ---------------------------------------------------------------------------
# Plot B: AUROC bar chart, models grouped by pathogen
# ---------------------------------------------------------------------------
bar_df = results_df.sort_values(["patho", "strain", "assay"]).reset_index(drop=True)
x = np.arange(len(bar_df))
colors = [MEASURE_COLORS[a] for a in bar_df["assay"]]
labels = [f"{p}\n{s}_{a}"
          for p, s, a in zip(bar_df["patho"], bar_df["strain"], bar_df["assay"])]

_, axs = st.create_figure(1, 1)
ax = axs.next()
ax.bar(x, bar_df["auroc"], color=colors)
ax.axhline(0.5, linestyle="--", color=nc.gray)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha="right")
ax.set_ylim(0, 1)
ax.legend(handles=[Patch(color=c, label=a) for a, c in MEASURE_COLORS.items()])
st.label(ax, xlabel="", ylabel="AUROC", title="EU-OPENSCREEN validation")

bar_path = os.path.join(output_dir, "12_auroc_barplot.png")
st.save_figure(bar_path)
print(f"Saved: {os.path.basename(bar_path)}")

print("\nDone.")
