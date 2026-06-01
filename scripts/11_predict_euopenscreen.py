"""
Predict the EU-OPENSCREEN library against CoADD models.

Constrained to:
  - inhib cutoff 50 only
  - mic cutoff 25 only

For each matching model:
  1. Predict all ~106k EU-OPENSCREEN SMILES → y_hat, y_score, y_rank
  2. Save flat JSON: {"smiles": [...], "y_hat": [...], "y_score": [...], "y_rank": [...]}
     (same format as 07 sample JSONs, extended with y_rank)
  3. Skip if prediction JSON already exists (checkpoint)

Per-pathogen plot: score distributions and pairwise rank scatter across models.

Output: output/11_predict_euopenscreen/
  11_{model_key}_predictions.json
  11_{patho}_model_comparison.png

Run with:
  conda run -n lazyqsar python scripts/11_predict_euopenscreen.py
"""

import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import stylia as st

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from lazyqsar.qsar import LazyClassifierQSAR

# Format: print | Style: ersilia
st.set_format("print")
st.set_style("ersilia")
nc = st.NamedColors()

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
EUOPENSCREEN_PATHOS = {
    "abaumannii", "calbicans", "ecoli", "efaecium",
    "kpneumoniae", "paeruginosa", "saureus",
}

CUTOFF_FILTER = {"inhib": "50", "mic": "25"}

euos_dir   = os.path.join(root, "..", "data", "raw", "eu-openscreen")
models_dir = os.path.join(root, "..", "output", "models")
output_dir = os.path.join(root, "..", "output", "11_predict_euopenscreen")
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Load EU-OPENSCREEN SMILES
# ---------------------------------------------------------------------------
all_smiles_df = pd.read_csv(os.path.join(euos_dir, "02_all_smiles.csv"))
all_smiles    = all_smiles_df["smiles"].tolist()
print(f"EU-OPENSCREEN library: {len(all_smiles):,} SMILES")

# ---------------------------------------------------------------------------
# Discover matching models (inhib_50 and mic_25 only)
# ---------------------------------------------------------------------------
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


matching = []
for model_key in sorted(os.listdir(models_dir)):
    try:
        patho, strain, assay, cutoff = parse_model_key(model_key)
    except (IndexError, ValueError):
        continue
    if patho not in EUOPENSCREEN_PATHOS:
        continue
    if cutoff != CUTOFF_FILTER.get(assay):
        continue
    model_path = os.path.join(models_dir, model_key, "model")
    if not os.path.exists(model_path) or not os.listdir(model_path):
        print(f"  [{model_key}] Model artefacts missing — skipping.")
        continue
    matching.append((patho, strain, assay, cutoff, model_key, model_path))

print(f"Matching models (inhib_50 / mic_25 only): {len(matching)}")
for entry in matching:
    print(entry)

# ---------------------------------------------------------------------------
# Predict and save flat JSON per model
# ---------------------------------------------------------------------------
for patho, strain, assay, cutoff, model_key, model_path in matching:
    pred_path = os.path.join(output_dir, f"11_{model_key}_predictions.json")

    if os.path.exists(pred_path):
        print(f"\n[{model_key}] Predictions already exist, skipping.")
        continue

    print(f"\n[{model_key}] Predicting {len(all_smiles):,} SMILES...")
    model = LazyClassifierQSAR.load(model_path)

    y_hat   = model.predict_proba(smiles_list=all_smiles)[:, 1].tolist()
    y_score = model.predict_score(smiles_list=all_smiles)[:, 1].tolist()
    y_rank  = model.predict_rank(smiles_list=all_smiles)[:, 1].tolist()

    result = {
        "smiles":  all_smiles,
        "y_hat":   y_hat,
        "y_score": y_score,
        "y_rank":  y_rank,
    }
    with open(pred_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"  Saved → {os.path.basename(pred_path)}")

# ---------------------------------------------------------------------------
# Per-pathogen comparison plots
# ---------------------------------------------------------------------------
patho_models = defaultdict(list)
for patho, strain, assay, cutoff, model_key, _ in matching:
    pred_path = os.path.join(output_dir, f"11_{model_key}_predictions.json")
    if os.path.exists(pred_path):
        patho_models[patho].append({
            "model_key": model_key,
            "assay":     assay,
            "pred_path": pred_path,
        })


def plot_score_histogram(ax, y_score, title, color):
    ax.hist(y_score, bins=50, color=color, edgecolor="none", alpha=0.8)
    st.label(ax, xlabel="y_score", ylabel="Count", title=title)


def plot_rank_scatter(ax, y_rank_a, y_rank_b, label_a, label_b):
    combined = np.array(y_rank_a) + np.array(y_rank_b)
    cm = st.FadingColormap("plum")
    cm.fit(combined)
    colors = cm.transform(combined)
    ax.scatter(y_rank_a, y_rank_b, c=colors, alpha=0.15, edgecolors="none",
               rasterized=True)
    st.label(ax, xlabel=label_a, ylabel=label_b, title="Rank comparison")


palette_colors = [nc.plum, nc.mint, nc.orange, nc.blue]

for patho, models in sorted(patho_models.items()):
    n = len(models)
    n_panels = n + (1 if n >= 2 else 0)
    _, axs = st.create_figure(1, n_panels, width=0.5)

    loaded = []
    for m in models:
        with open(m["pred_path"]) as fh:
            data = json.load(fh)
        loaded.append((m, data))

    for i, (m, data) in enumerate(loaded):
        ax = axs.next()
        plot_score_histogram(ax, data["y_score"], m["model_key"],
                             palette_colors[i % len(palette_colors)])

    if n >= 2:
        ax = axs.next()
        m_a, data_a = loaded[0]
        m_b, data_b = loaded[1]
        plot_rank_scatter(
            ax,
            data_a["y_rank"],
            data_b["y_rank"],
            label_a=m_a["model_key"],
            label_b=m_b["model_key"],
        )

    out = os.path.join(output_dir, f"11_{patho}_model_comparison.png")
    st.save_figure(out)
    print(f"Saved: {os.path.basename(out)}")

print("\nDone.")
