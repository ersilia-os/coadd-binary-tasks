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

Per-pathogen plot: one score histogram per strain (inhib/mic measures overlaid)
and, below each, an inhib-vs-mic rank scatter for that same strain.

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

# Scatter plots overplot badly with ~106k points; uniformly subsample for display.
SCATTER_SAMPLE = 10000
SCATTER_SEED   = 42

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
            "strain":    strain,
            "assay":     assay,
            "cutoff":    cutoff,
            "pred_path": pred_path,
        })

# One colour per measure, consistent across strains.
MEASURE_COLORS = {"inhib": nc.blue, "mic": nc.orange}


def plot_score_histograms(ax, measures, title):
    """Overlay the y_score distribution of each measure (inhib/mic) for a strain."""
    for label, data, color in measures:
        ax.hist(data["y_score"], bins=50, color=color, edgecolor="none",
                alpha=0.6, label=label)
    ax.legend()
    st.label(ax, xlabel="y_score", ylabel="Count", title=title)


def plot_rank_scatter(ax, y_rank_x, y_rank_y, label_x, label_y):
    x = np.array(y_rank_x)
    y = np.array(y_rank_y)

    # Uniformly subsample for display — the same indices on both axes.
    if x.size > SCATTER_SAMPLE:
        rng = np.random.default_rng(SCATTER_SEED)
        idx = rng.choice(x.size, size=SCATTER_SAMPLE, replace=False)
        x, y = x[idx], y[idx]

    ax.scatter(x, y, color=nc.plum, s=st.MARKERSIZE_SMALL, alpha=0.4,
               edgecolors="none", rasterized=True)
    st.label(ax, xlabel=label_x, ylabel=label_y)


for patho, models in sorted(patho_models.items()):
    # Group models by strain, preserving discovery order.
    strains = []
    by_strain = defaultdict(dict)   # strain -> {assay: model}
    for m in models:
        if m["strain"] not in by_strain:
            strains.append(m["strain"])
        by_strain[m["strain"]][m["assay"]] = m

    cache = {}
    for m in models:
        with open(m["pred_path"]) as fh:
            cache[m["model_key"]] = json.load(fh)

    # Each strain gets a histogram (its measures overlaid). A strain also gets a
    # rank scatter — inhib vs mic, directly below its histogram — only when both
    # measures exist. We never compare across strains.
    has_scatter = any({"inhib", "mic"} <= set(by_strain[s]) for s in strains)
    nrows = 2 if has_scatter else 1
    # A single strain is one narrow column — don't stretch it to full width.
    fig_kw = {"width": 0.5} if len(strains) == 1 else {}
    _, axs = st.create_figure(nrows, len(strains), **fig_kw)

    # Row 1: one histogram per strain, measures overlaid with a legend.
    for s in strains:
        ax = axs.next()
        measures = []
        for assay in ("inhib", "mic"):
            if assay in by_strain[s]:
                m = by_strain[s][assay]
                measures.append((f"{assay}_{m['cutoff']}",
                                 cache[m["model_key"]], MEASURE_COLORS[assay]))
        plot_score_histograms(ax, measures, title=s)

    # Row 2: one inhib-vs-mic scatter per strain (blank where a measure is missing).
    if has_scatter:
        for s in strains:
            ax = axs.next()
            if {"inhib", "mic"} <= set(by_strain[s]):
                m_i = by_strain[s]["inhib"]
                m_m = by_strain[s]["mic"]
                plot_rank_scatter(
                    ax,
                    cache[m_i["model_key"]]["y_rank"],
                    cache[m_m["model_key"]]["y_rank"],
                    label_x=f"inhib_{m_i['cutoff']} rank",
                    label_y=f"mic_{m_m['cutoff']} rank",
                )
            else:
                ax.axis("off")

    out = os.path.join(output_dir, f"11_{patho}_model_comparison.png")
    st.save_figure(out)
    print(f"Saved: {os.path.basename(out)}")

print("\nDone.")
