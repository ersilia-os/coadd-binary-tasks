"""
Train LazyQSAR models for CoADD cytotoxicity datasets and compare predictions
against bioactivity models (script 07) on a random ChEMBL sample.

Part 1 — Train cytotoxicity models
  CC50 and HC10 at cutoffs [10, 25, 50] µM.
  5-fold cross-validation + final model; skips cutoffs with < MIN_POSITIVES actives.
  Outputs:
    output/models/cc50_{cutoff}/crossval_report.json  (+ model/)
    output/models/hc10_{cutoff}/crossval_report.json  (+ model/)
    output/08_lq_cytotox/cv_summary.csv
    output/08_lq_cytotox/08_cc50_performance.png
    output/08_lq_cytotox/08_hc10_performance.png

Part 2 — Predict on random ChEMBL sample
  Each trained cytotoxicity model predicts on data/raw/chembl_random.csv.
  Outputs:
    output/08_lq_cytotox/08_sample_{model}.json  (per model)
    output/08_lq_cytotox/08_cytotox_predictions.csv

Part 3 — Spearman correlation: cytotoxicity vs bioactivity
  Cross-correlate cytotoxicity and bioactivity predictions (from
  07_model_prediction_correlation.csv) on the same ChEMBL sample.
  Outputs:
    output/08_lq_cytotox/08_cytotox_bioactivity_correlation.csv
    output/08_lq_cytotox/08_cytotox_bioactivity_correlation.png

Run with:
  conda run -n lazyqsar python scripts/08_lq_cytotox.py

Parallel runs: train subsets with --only KEY[,KEY...] --train-only (KEY e.g. cc50_25),
then run once without flags for the summary, plots and predictions. Part 3 reads
step 07's 07_model_prediction_correlation.csv, so the final run goes after step 07.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import stylia as st
from scipy.stats import spearmanr

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from lazyqsar.qsar import LazyClassifierQSAR
from src.lazyqsar_utils import run_crossval, predict_and_save
from src.plotting_utils import (
    plot_class_balance,
    plot_cross_corr,
    plot_roc_folds,
    plot_scores,
)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
N_FOLDS       = 5
MIN_POSITIVES = 20
CUTOFFS       = [10, 25, 50]

parser = argparse.ArgumentParser(description="Train LazyQSAR cytotoxicity models.")
parser.add_argument("--only", default=None,
                    help="Comma-separated model keys to train; all others are skipped")
parser.add_argument("--train-only", action="store_true",
                    help="Stop after training, before the summary, plots and predictions")
args = parser.parse_args()
ONLY = set(args.only.split(",")) if args.only else None

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
cytotox_dir = os.path.join(root, "..", "data", "processed", "coadd", "06_cytotox")
models_dir  = os.path.join(root, "..", "output", "models")
output_dir  = os.path.join(root, "..", "output", "08_lq_cytotox")
bioact_csv  = os.path.join(
    root, "..", "output", "07_lq_coadd", "07_model_prediction_correlation.csv"
)

os.makedirs(models_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

DATASETS = {
    "cc50": {"file": "cc50.csv", "col_prefix": "cc50"},
    "hc10": {"file": "hc10.csv", "col_prefix": "hc10"},
}

# ---------------------------------------------------------------------------
# Part 1 — Train cytotoxicity models
# ---------------------------------------------------------------------------
print("=" * 60)
print("PART 1 — Train cytotoxicity models")
print("=" * 60)

summary_rows    = []
available_models = []  # list of {key, assay, cutoff, report_path}

for assay_name, cfg in DATASETS.items():
    df     = pd.read_csv(os.path.join(cytotox_dir, cfg["file"]))
    prefix = cfg["col_prefix"]

    print(f"\n{'='*60}")
    print(f"{assay_name.upper()}  ({len(df):,} compounds)")
    print(f"{'='*60}")

    for cutoff in CUTOFFS:
        col = f"{prefix}_{cutoff}"
        if col not in df.columns:
            print(f"  [{col}] Column missing, skipping.")
            continue

        df_def = df[df[col] != -1][["std_smiles", col]].dropna()
        smiles = df_def["std_smiles"].tolist()
        y      = df_def[col].astype(int).tolist()

        n_tot = len(y)
        n_pos = int(sum(y))
        n_neg = n_tot - n_pos

        if ONLY is not None and f"{assay_name}_{cutoff}" not in ONLY:
            continue

        if n_pos < MIN_POSITIVES:
            print(f"  [{col}] Too few positives ({n_pos}), skipping.")
            continue

        print(f"\n  [{col}]  n={n_tot:,}  positives={n_pos} ({n_pos/n_tot:.1%})")

        model_key = f"{assay_name}_{cutoff}"
        model_dir = os.path.join(models_dir, model_key)
        os.makedirs(model_dir, exist_ok=True)

        # Cross-validation
        cv_path = os.path.join(model_dir, "crossval_report.json")
        if os.path.exists(cv_path):
            print(f"  [{col}] CV report exists, loading.")
            with open(cv_path) as fh:
                report = json.load(fh)
        else:
            report = run_crossval(
                smiles_list=smiles,
                y=y,
                n_folds=N_FOLDS,
                output_path=cv_path,
                mode="slow",
            )
        fold_aucs = [v["roc_auc"] for v in report.values()]
        mean_auc  = float(np.mean(fold_aucs))
        std_auc   = float(np.std(fold_aucs))
        print(f"  CV AUC: {mean_auc:.4f} ± {std_auc:.4f}")

        # Final model
        final_model_dir = os.path.join(model_dir, "model")
        if os.path.exists(final_model_dir) and os.listdir(final_model_dir):
            print(f"  [{col}] Final model exists, skipping.")
        else:
            os.makedirs(final_model_dir, exist_ok=True)
            model = LazyClassifierQSAR(mode="slow")
            model.fit(smiles_list=smiles, y=np.array(y))
            model.save(final_model_dir)
            print(f"  Final model saved → output/models/{model_key}/model/")

        summary_rows.append({
            "assay":       assay_name,
            "cutoff":      cutoff,
            "column":      col,
            "n_compounds": n_tot,
            "n_active":    n_pos,
            "n_inactive":  n_neg,
            "active_rate": round(n_pos / n_tot, 4),
            "cv_auc_mean": round(mean_auc, 4),
            "cv_auc_std":  round(std_auc, 4),
        })
        available_models.append({
            "key":         model_key,
            "assay":       assay_name,
            "cutoff":      cutoff,
            "report_path": cv_path,
        })

if args.train_only:
    print("\nTraining done (--train-only).")
    sys.exit(0)

if summary_rows:
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(output_dir, "cv_summary.csv"), index=False
    )
    print(f"\nSaved: output/08_lq_cytotox/cv_summary.csv")

# Performance plots — one figure per assay (rows = cutoffs, cols = balance/ROC/scores)
for assay_name in DATASETS:
    assay_models = [m for m in available_models if m["assay"] == assay_name]
    if not assay_models:
        continue
    n = len(assay_models)
    _, axs = st.create_figure(n, 3, width_ratios=[0.3, 1, 1])
    for m in assay_models:
        with open(m["report_path"]) as fh:
            data = json.load(fh)
        plot_class_balance(axs.next(), data, title=m["key"])
        ax = axs.next()
        plot_roc_folds(ax, data)
        st.label(ax, xlabel="FPR", ylabel="TPR", title=m["key"])
        plot_scores(axs.next(), data, title=m["key"])
    out = os.path.join(output_dir, f"08_{assay_name}_performance.png")
    st.save_figure(out)
    print(f"Saved: {os.path.basename(out)}")

# ---------------------------------------------------------------------------
# Part 2 — Predict on random ChEMBL sample
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("PART 2 — Predict on random ChEMBL sample")
print("=" * 60)

sample_smiles = (
    pd.read_csv(os.path.join(root, "..", "data", "raw", "chembl_random.csv"))
    ["smiles"].dropna().tolist()
)
print(f"ChEMBL sample: {len(sample_smiles):,} SMILES")

cytotox_PRED_CSV = os.path.join(output_dir, "08_cytotox_predictions.csv")

if os.path.exists(cytotox_PRED_CSV):
    print(f"Loading cached predictions from {os.path.basename(cytotox_PRED_CSV)}")
    cytotox_pred_df = pd.read_csv(cytotox_PRED_CSV, index_col=0)
else:
    pred_dict = {}
    for m in available_models:
        final_model_dir = os.path.join(models_dir, m["key"], "model")
        if not os.path.exists(final_model_dir):
            print(f"  [{m['key']}] Final model not found — skipping.")
            continue
        print(f"  Predicting with {m['key']}...")
        model = LazyClassifierQSAR.load(final_model_dir)
        pred_path = os.path.join(output_dir, f"08_sample_{m['key']}.json")
        pred_dict[m["key"]] = predict_and_save(model, sample_smiles, pred_path)

    cytotox_pred_df = pd.DataFrame(pred_dict, index=sample_smiles)
    cytotox_pred_df.to_csv(cytotox_PRED_CSV)
    print(f"Saved: {os.path.basename(cytotox_PRED_CSV)}")

# ---------------------------------------------------------------------------
# Part 3 — Spearman correlation: cytotoxicity vs bioactivity
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("PART 3 — Spearman correlation: cytotoxicity vs bioactivity")
print("=" * 60)

bioact_pred_df = pd.read_csv(bioact_csv, index_col=0)
print(f"Bioactivity models: {len(bioact_pred_df.columns)}")
print(f"cytotox models:     {len(cytotox_pred_df.columns)}")

shared_idx = cytotox_pred_df.index.intersection(bioact_pred_df.index)
print(f"Shared compounds:   {len(shared_idx):,}")

cytotox_aligned = cytotox_pred_df.loc[shared_idx]
bioact_aligned  = bioact_pred_df.loc[shared_idx]

# Cross-correlation: rows = cytotox models, cols = bioactivity models
corr_rows = {}
for ct_col in cytotox_aligned.columns:
    row = {}
    for ba_col in bioact_aligned.columns:
        r, _ = spearmanr(cytotox_aligned[ct_col], bioact_aligned[ba_col])
        row[ba_col] = round(float(r), 4)
    corr_rows[ct_col] = row

corr_df = pd.DataFrame(corr_rows).T  # (n_cytotox × n_bioact)
corr_df.to_csv(os.path.join(output_dir, "08_cytotox_bioactivity_correlation.csv"))
print(f"\nSpearman r (cytotox vs bioactivity):")
print(corr_df.to_string())
print(f"\nSaved: 08_cytotox_bioactivity_correlation.csv")

# Heatmap
_, axs = st.create_figure(1, 1)
plot_cross_corr(
    axs.next(),
    corr_df,
    title="Spearman r — cytotoxicity vs bioactivity",
)
out = os.path.join(output_dir, "08_cytotox_bioactivity_correlation.png")
st.save_figure(out)
print(f"Saved: {os.path.basename(out)}")

print("\nDone.")
