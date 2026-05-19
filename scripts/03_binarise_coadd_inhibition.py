"""
Binarise CoADD inhibition data per pathogen and strain.

For each pathogen in COADD_INH_PATHO:
  1. Load preprocessed file from 02_preprocess_inh/
  2. Exclude measurements from 300 ug/mL (≈ 832.2 µM)
  3. For duplicate SMILES within the same scope: average the raw inhibition
     value across all retained measurements, then binarise the average
  4. Produce per-strain files and one merged (all-strains) file per pathogen
  5. Skip any output file with < MIN_COMPOUNDS unique SMILES

Outputs:
  data/processed/coadd/03_binarised_inhibition/{patho_code}_{strain_code}.csv
  data/processed/coadd/03_binarised_inhibition/{patho_code}_merged.csv
  Columns: std_smiles, inchikey, mw, inhib_50, inhib_75, inhib_90

  output/03_binarise_coadd_inhibition/summary.csv
  output/03_binarise_coadd_inhibition/{patho_code}_binarised.png
"""

import os
import sys
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.default import COADD_INH_PATHO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
preprocess_dir = os.path.join(root, "..", "data", "processed", "coadd", "02_preprocess_inh")
binarised_dir  = os.path.join(root, "..", "data", "processed", "coadd", "03_binarised_inhibition")
output_dir     = os.path.join(root, "..", "output", "03_binarise_coadd_inhibition")
config_dir     = os.path.join(root, "..", "config", "manual")

os.makedirs(binarised_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Cutoffs from config
# ---------------------------------------------------------------------------
cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))
inhib_cuts = cutoffs_cfg[cutoffs_cfg["assay_type"] == "inhib"]
CUTOFFS = sorted(set(
    inhib_cuts[["cutoff_low", "cutoff_mid", "cutoff_high"]].values.flatten().tolist()
))
CUTOFF_COLS = [f"inhib_{c:.0f}" for c in CUTOFFS]
print(f"Cutoffs: {CUTOFFS}  →  columns: {CUTOFF_COLS}")

# ---------------------------------------------------------------------------
# 300 ug/mL exclusion: round concentration_um to 1 decimal, exclude 832.2 µM
# ---------------------------------------------------------------------------
EXCLUDE_CONC_UM = 832.2


def drop_300ugml(df):
    return df[df["concentration_um"].round(1) != EXCLUDE_CONC_UM].copy()


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------
def binarise(df):
    """Average value per std_smiles, binarise at each cutoff.

    Returns a DataFrame with columns:
        std_smiles, inchikey, mw, inhib_50, inhib_75, inhib_90
    Only rows with operator == "=" are used in the average.
    """
    df = df[df["operator"] == "="].copy()
    if df.empty:
        return pd.DataFrame(columns=["std_smiles", "inchikey", "mw"] + CUTOFF_COLS)

    # Average inhibition value per unique SMILES, with replica count and std
    agg = (
        df.groupby("std_smiles")["value"]
        .agg(avg_value="mean", replicas="count", std="std")
        .reset_index()
    )
    agg["std"] = agg["std"].fillna(0)  # std is NaN when replicas == 1
    avg_vals = agg

    # Carry forward inchikey, mw, and aggregated operator (first occurrence per SMILES)
    meta = (
        df[["std_smiles", "inchikey", "mw"]]
        .drop_duplicates(subset="std_smiles")
        .reset_index(drop=True)
    )
    op_agg = (
        df.groupby("std_smiles")["operator"]
        .apply(lambda ops: ops.iloc[0] if ops.nunique() == 1 else "mixed")
        .reset_index()
    )

    out = avg_vals.merge(meta, on="std_smiles", how="left").merge(op_agg, on="std_smiles", how="left")

    for cutoff, col in zip(CUTOFFS, CUTOFF_COLS):
        out[col] = (out["avg_value"] >= cutoff).astype(int)

    out["avg_value"] = out["avg_value"].round(3)
    out["std"]       = out["std"].round(3)
    out["mw"]        = out["mw"].round(3)

    return out[["std_smiles", "inchikey", "mw", "avg_value", "std", "operator", "replicas"] + CUTOFF_COLS].rename(
        columns={"avg_value": "value"}
    ).reset_index(drop=True)


def save_dataset(df_bin, filepath, label):
    """Save dataset and return summary row."""
    n = len(df_bin)
    df_bin.to_csv(filepath, index=False)
    row = {"label": label, "n_compounds": n}
    for col in CUTOFF_COLS:
        n_active = int(df_bin[col].sum())
        row[f"n_active_{col}"] = n_active
        row[f"active_rate_{col}"] = round(n_active / n, 4)
    print(f"  SAVED {label}: {n} compounds | "
          + " | ".join(f"{col}={row[f'active_rate_{col}']:.1%}" for col in CUTOFF_COLS))
    return row


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
summary_rows = []

for patho_code in COADD_INH_PATHO:
    fpath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    if not os.path.exists(fpath):
        print(f"\n[{patho_code}] Preprocessed file not found, skipping.")
        continue

    df = pd.read_csv(fpath)
    print(f"\n[{patho_code}] {len(df):,} rows loaded")

    # Step 1: exclude 300 ug/mL for all downstream processing
    df_filtered = drop_300ugml(df)
    print(f"  {len(df_filtered):,} rows after excluding 300 ug/mL ({len(df) - len(df_filtered):,} dropped)")

    patho_summary = []

    # Step 2: binarise per strain, averaging duplicates within each strain
    for strain_code in sorted(df_filtered["strain_code"].dropna().unique()):
        df_s = df_filtered[df_filtered["strain_code"] == strain_code]
        df_bin = binarise(df_s)
        fname = f"{patho_code}_{strain_code}.csv"
        row = save_dataset(df_bin, os.path.join(binarised_dir, fname),
                           f"{patho_code}/{strain_code}")
        row["patho_code"] = patho_code
        row["strain_code"] = strain_code
        row["file"] = fname
        summary_rows.append(row)
        patho_summary.append({"label": strain_code, "n": len(df_bin),
                               **{col: df_bin[col].mean() for col in CUTOFF_COLS}})

    # Step 3: merged — only when more than one strain is present
    strains = df_filtered["strain_code"].dropna().unique()
    if len(strains) > 1:
        df_bin_merged = binarise(df_filtered)
        fname = f"{patho_code}_merged.csv"
        row = save_dataset(df_bin_merged, os.path.join(binarised_dir, fname),
                           f"{patho_code}/merged")
        row["patho_code"] = patho_code
        row["strain_code"] = "merged"
        row["file"] = fname
        summary_rows.append(row)
        patho_summary.append({"label": "merged", "n": len(df_bin_merged),
                               **{col: df_bin_merged[col].mean() for col in CUTOFF_COLS}})
    else:
        print(f"  SKIP merged: only 1 strain present")

    # -----------------------------------------------------------------------
    # Plot: dataset size (A) + active fraction by cutoff (B)
    # -----------------------------------------------------------------------
    if not patho_summary:
        continue

    import stylia
    from matplotlib import pyplot as plt

    # Format: print | Style: ersilia
    stylia.set_format("print")
    stylia.set_style("ersilia")

    labels   = [r["label"] for r in patho_summary]
    ns       = [r["n"] for r in patho_summary]
    rates    = {col: [r[col] for r in patho_summary] for col in CUTOFF_COLS}
    x        = np.arange(len(labels))
    n_groups = len(CUTOFFS)
    bar_w    = 0.25

    pal    = stylia.CategoricalPalette("ersilia")
    colors = pal.get(n_groups)
    nc     = stylia.NamedColors()

    fig_width = 0.5 if len(labels) == 1 else 1.0
    fig, axs = stylia.create_figure(1, 2, width=fig_width)

    # Panel A — compound count (vertical bar)
    ax = axs.next()
    ax.bar(x, ns, color=nc.blue)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    stylia.label(ax, ylabel="Compounds", xlabel="", title=patho_code)

    # Panel B — active fraction per cutoff (vertical grouped bar)
    ax = axs.next()
    for i, (col, cutoff) in enumerate(zip(CUTOFF_COLS, CUTOFFS)):
        offset = (i - (n_groups - 1) / 2) * bar_w
        ax.bar(x + offset, rates[col], width=bar_w, color=colors[i],
               label=f"≥{cutoff:.0f}%")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.legend()
    stylia.label(ax, ylabel="Active fraction", xlabel="", title=patho_code)

    stylia.save_figure(os.path.join(output_dir, f"{patho_code}_binarised.png"))
    print(f"  → plot saved: {patho_code}_binarised.png")

# ---------------------------------------------------------------------------
# Summary CSV
# ---------------------------------------------------------------------------
if summary_rows:
    summary_df = pd.DataFrame(summary_rows)[
        ["patho_code", "strain_code", "n_compounds"]
        + [f"n_active_{col}" for col in CUTOFF_COLS]
        + [f"active_rate_{col}" for col in CUTOFF_COLS]
    ]
    summary_df.to_csv(os.path.join(output_dir, "summary.csv"), index=False)
    print(f"\nSummary written: {len(summary_df)} rows → output/03_binarise_coadd_inhibition/summary.csv")

print("\nDone.")
