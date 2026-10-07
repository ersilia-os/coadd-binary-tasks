"""
Train LazyQSAR models for all selected CoADD datasets.

For each dataset with keep=True in the selection_table (output/06_coadd_analysis/),
and for each activity cutoff (low / mid / high):
  1. Load the binarised file from data/processed/coadd/
  2. Drop inconclusive labels (-1); skip if fewer than MIN_POSITIVES active compounds
  3. Run 5-fold cross-validation (StratifiedShuffleSplit) → crossval_report.json
  4. Train final model on all data (mode="slow") → saved with model.save()

Output structure:
  output/models/{patho}_{strain}_{assay}_{cutoff_value}/
    crossval_report.json   — per-fold y_true, y_hat, y_score, y_rank, roc_auc
    model/                 — final LazyQSAR model artefacts (from model.save())

  output/07_lq_coadd/cv_summary.csv  — one row per (dataset, cutoff) with mean CV AUC

Run with:
  conda run -n lazyqsar python scripts/07_lq_coadd.py

Parallel runs: train disjoint subsets with --only KEY[,KEY...] --train-only (KEY is the
model folder name, e.g. saureus_ATCC43300_inhib_50), then run once without flags to
collect every CV report into the summary, plots and ChEMBL-sample predictions.
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
from collections import defaultdict
import stylia as st

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from lazyqsar.qsar import LazyClassifierQSAR
from src.lazyqsar_utils import run_crossval, predict_and_save
from src.plotting_utils import plot_class_balance, plot_roc_folds, plot_scores, plot_corr_matrix, plot_overlap_matrix, compute_hit_overlap

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
N_FOLDS       = 5
MIN_POSITIVES = 20   # skip a cutoff if fewer active compounds than this

parser = argparse.ArgumentParser(description="Train LazyQSAR models for CoADD datasets.")
parser.add_argument("--only", default=None,
                    help="Comma-separated model keys to train; all others are skipped")
parser.add_argument("--train-only", action="store_true",
                    help="Stop after training, before the summary, plots and predictions")
args = parser.parse_args()
ONLY = set(args.only.split(",")) if args.only else None

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
config_dir  = os.path.join(root, "..", "config", "manual")
inh_bin_dir = os.path.join(root, "..", "data", "processed", "coadd", "03_binarised_inhibition")
mic_bin_dir = os.path.join(root, "..", "data", "processed", "coadd", "05_binarised_mic")
models_dir  = os.path.join(root, "..", "output", "models")
output_dir  = os.path.join(root, "..", "output", "07_lq_coadd")

os.makedirs(models_dir, exist_ok=True)
os.makedirs(output_dir,  exist_ok=True)

# ---------------------------------------------------------------------------
# Load selection table (kept datasets only)
# ---------------------------------------------------------------------------
sel = pd.read_csv(
    os.path.join(root, "..", "output", "06_coadd_analysis", "selection_table.csv")
)
sel = sel[sel["keep"]].reset_index(drop=True)

print(f"Datasets to train: {len(sel)}")

# ---------------------------------------------------------------------------
# Build cutoff tier → (column_name, cutoff_value) map per assay type
# ---------------------------------------------------------------------------
cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))


def get_cutoff_map(assay_type):
    """Returns dict: tier → (df_column, cutoff_value)."""
    row    = cutoffs_cfg[cutoffs_cfg["assay_type"] == assay_type].iloc[0]
    prefix = "inhib" if assay_type == "inhib" else "mic"
    return {
        "low":  (f"{prefix}_{int(row['cutoff_low'])}",  int(row["cutoff_low"])),
        "mid":  (f"{prefix}_{int(row['cutoff_mid'])}",  int(row["cutoff_mid"])),
        "high": (f"{prefix}_{int(row['cutoff_high'])}", int(row["cutoff_high"])),
    }


# ---------------------------------------------------------------------------
# Main training loop
# ---------------------------------------------------------------------------
summary_rows = []

for _, ds in sel.iterrows():
    patho  = ds["patho_code"]
    strain = ds["strain_code"]
    assay  = ds["assay_type"]

    bin_dir = inh_bin_dir if assay == "inhib" else mic_bin_dir
    fpath   = os.path.join(bin_dir, f"{patho}_{strain}.csv")
    if not os.path.exists(fpath):
        print(f"\n[{patho}/{strain}/{assay}] Binarised file not found, skipping.")
        continue

    df = pd.read_csv(fpath)
    cutoff_map = get_cutoff_map(assay)

    print(f"\n{'='*60}")
    print(f"{patho} / {strain} / {assay}  ({len(df):,} compounds)")
    print(f"{'='*60}")

    for tier, (col, cutoff_val) in cutoff_map.items():
        if col not in df.columns:
            print(f"  [{tier}] Column {col} missing, skipping.")
            continue

        # Drop inconclusive (-1); keep only 0/1
        df_def = df[df[col] != -1][["std_smiles", col]].dropna()
        smiles = df_def["std_smiles"].tolist()
        y      = df_def[col].astype(int).tolist()

        n_tot = len(y)
        n_pos = int(sum(y))
        n_neg = n_tot - n_pos

        if ONLY is not None and f"{patho}_{strain}_{assay}_{cutoff_val}" not in ONLY:
            continue

        if n_pos < MIN_POSITIVES:
            print(f"  [{tier} / {col}] Too few positives ({n_pos}), skipping.")
            continue

        print(f"\n  [{tier} / {col}]  n={n_tot:,}  positives={n_pos} ({n_pos/n_tot:.1%})")

        model_dir = os.path.join(
            models_dir, f"{patho}_{strain}_{assay}_{cutoff_val}"
        )
        os.makedirs(model_dir, exist_ok=True)

        # --- Cross-validation ---
        cv_path = os.path.join(model_dir, "crossval_report.json")
        if os.path.exists(cv_path):
            print(f"  [{tier} / {col}] CV report already exists, loading.")
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

        # --- Final model (all data) ---
        final_model_dir = os.path.join(model_dir, "model")
        if os.path.exists(final_model_dir) and os.listdir(final_model_dir):
            print(f"  [{tier} / {col}] Final model already exists, skipping.")
        else:
            os.makedirs(final_model_dir, exist_ok=True)
            model = LazyClassifierQSAR(mode="slow")
            model.fit(smiles_list=smiles, y=np.array(y))
            model.save(final_model_dir)
            print(f"  Final model saved → output/models/{patho}_{strain}_{assay}_{cutoff_val}/model/")

        summary_rows.append({
            "patho_code":   patho,
            "strain_code":  strain,
            "assay_type":   assay,
            "cutoff_tier":  tier,
            "cutoff_value": cutoff_val,
            "column":       col,
            "n_compounds":  n_tot,
            "n_active":     n_pos,
            "n_inactive":   n_neg,
            "active_rate":  round(n_pos / n_tot, 4),
            "cv_auc_mean":  round(mean_auc, 4),
            "cv_auc_std":   round(std_auc, 4),
        })

if args.train_only:
    print("\nTraining done (--train-only).")
    sys.exit(0)

# ---------------------------------------------------------------------------
# Summary CSV
# ---------------------------------------------------------------------------
if summary_rows:
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(output_dir, "cv_summary.csv"), index=False)
    print(f"\n{'='*60}")
    print("CV SUMMARY")
    print(f"{'='*60}")
    print(summary_df[["patho_code", "strain_code", "assay_type",
                       "cutoff_tier", "n_compounds", "n_active",
                       "cv_auc_mean", "cv_auc_std"]].to_string(index=False))
    print(f"\nSaved: output/07_lq_coadd/cv_summary.csv")

print("\nDone.")

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Collect all datasets that have a crossval report
available_datasets = []
for _, ds in sel.iterrows():
    patho  = ds["patho_code"]
    strain = ds["strain_code"]
    assay  = ds["assay_type"]
    for tier, (col, cutoff_val) in get_cutoff_map(assay).items():
        dataset     = f"{patho}_{strain}_{assay}_{cutoff_val}"
        report_path = os.path.join(models_dir, dataset, "crossval_report.json")
        if os.path.exists(report_path):
            available_datasets.append({
                "dataset":     dataset,
                "patho":       patho,
                "strain":      strain,
                "assay":       assay,
                "cutoff":      cutoff_val,
                "tier":        tier,
                "report_path": report_path,
            })

# Group by (patho, strain, assay); one figure per group (rows = cutoff tiers, max 3)
psa_groups = defaultdict(list)
for d in available_datasets:
    key = (d["patho"], d["strain"], d["assay"])
    psa_groups[key].append(d)

for (patho, strain, assay), datasets in psa_groups.items():
    n = len(datasets)
    _, axs = st.create_figure(n, 3, width_ratios=[0.3, 1, 1])
    for d in datasets:
        with open(d["report_path"]) as fh:
            data = json.load(fh)
        ax = axs.next()
        plot_class_balance(ax, data, title=d["dataset"])
        ax = axs.next()
        plot_roc_folds(ax, data)
        st.label(ax, xlabel="FPR", ylabel="TPR", title=d["dataset"])
        ax = axs.next()
        plot_scores(ax, data, title=d["dataset"])
    out = os.path.join(output_dir, f"07_{patho}_{strain}_{assay}_performance.png")
    st.save_figure(out)
    print(f"Saved: {out}")

# ---------------------------------------------------------------------------
# Cross-model prediction correlation on ChEMBL sample
# ---------------------------------------------------------------------------
CORR_CSV = os.path.join(output_dir, "07_model_prediction_correlation.csv")

sample_smiles = pd.read_csv(
    os.path.join(root, "..", "data", "raw", "chembl_random.csv")
)["smiles"].dropna().tolist()

if os.path.exists(CORR_CSV):
    print(f"Loading cached predictions from {os.path.basename(CORR_CSV)}")
    pred_df = pd.read_csv(CORR_CSV, index_col=0)
else:
    pred_dict = {}
    for d in available_datasets:
        key       = d["dataset"]
        model_dir = os.path.join(models_dir, key, "model")
        if not os.path.exists(model_dir):
            print(f"  [{key}] Final model not found — skipping.")
            continue
        print(f"  Predicting with {key}...")
        model    = LazyClassifierQSAR.load(model_dir)
        pred_path = os.path.join(output_dir, f"07_sample_{key}.json")
        pred_dict[key] = predict_and_save(model, sample_smiles, pred_path)

    if len(pred_dict) < 2:
        print("  Fewer than 2 models available — skipping correlation plots.")
    else:
        pred_df = pd.DataFrame(pred_dict, index=sample_smiles)
        pred_df.to_csv(CORR_CSV)

if "pred_df" in dir() and len(pred_df.columns) >= 2:
    _, axs = st.create_figure(1, 3)
    ax = axs.next()
    plot_corr_matrix(ax, pred_df.corr(method="spearman"),
                        title=f"Spearman r — all ({len(pred_df):,})", show_labels=False)
    ax = axs.next()
    plot_overlap_matrix(ax, compute_hit_overlap(pred_df, 1000),
                        title="Hit overlap — top 1,000", show_labels=False)
    ax = axs.next()
    plot_overlap_matrix(ax, compute_hit_overlap(pred_df, 100),
                        title="Hit overlap — top 100", show_labels=False)

    out = os.path.join(output_dir, "07_model_prediction_correlation.png")
    st.save_figure(out)
    print(f"Saved: {os.path.basename(out)}")
