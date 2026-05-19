"""
Binarise CoADD MIC data per pathogen and strain.

For each pathogen in COADD_MIC_PATHO:
  1. Load preprocessed file from 04_preprocess_mic/
  2. Binarise per strain, aggregating replicates
  3. Create merged (all-strains) file if more than one strain present
  4. Save plots and summary CSV

Binarisation logic (direction = -1, lower MIC = more active):
  operator "=" : 1 if value <= cutoff  else  0
  operator ">" : 0 if value >= cutoff  else  -1  (MIC > value; if value >= cutoff, definitely inactive)
  operator "<" : 1 if value <= cutoff  else  -1  (MIC < value; if value <= cutoff, definitely active)

Aggregation across replicates:
  - Compute value (mean of all numerics), std (0 if n=1), replicas (count)
  - For each cutoff, apply binarize_mic to each row → (1, 0, -1)
  - Use definitive results (0s and 1s) if any exist; threshold mean at 0.5
  - If all inconclusive (-1): final label = -1

Outputs:
  data/processed/coadd/05_binarised_mic/{patho_code}_{strain_code}.csv
  data/processed/coadd/05_binarised_mic/{patho_code}_merged.csv  (only if >1 strain)
  Columns: std_smiles, inchikey, mw, value, std, replicas, mic_10, mic_25, mic_50

  output/05_binarise_coadd_mic/summary.csv
  output/05_binarise_coadd_mic/{patho_code}_binarised.png
"""

import os
import sys
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.default import COADD_MIC_PATHO

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
preprocess_dir = os.path.join(root, "..", "data", "processed", "coadd", "04_preprocess_mic")
binarised_dir  = os.path.join(root, "..", "data", "processed", "coadd", "05_binarised_mic")
output_dir     = os.path.join(root, "..", "output", "05_binarise_coadd_mic")
config_dir     = os.path.join(root, "..", "config", "manual")

os.makedirs(binarised_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Cutoffs from config (MIC, direction = -1)
# ---------------------------------------------------------------------------
cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))
mic_cuts = cutoffs_cfg[cutoffs_cfg["assay_type"] == "mic"]
CUTOFFS = sorted(set(
    mic_cuts[["cutoff_low", "cutoff_mid", "cutoff_high"]].values.flatten().tolist()
))
CUTOFF_COLS = [f"mic_{c:.0f}" for c in CUTOFFS]
print(f"Cutoffs: {CUTOFFS}  →  columns: {CUTOFF_COLS}")


# ---------------------------------------------------------------------------
# Binarisation logic for MIC (direction = -1)
# ---------------------------------------------------------------------------
def binarize_mic(operator, value, cutoff):
    """Return 1 (active), 0 (inactive), or -1 (inconclusive).

    Lower MIC = more active (direction = -1).
      "=" : definitive — active if value <= cutoff
      ">" : MIC > value; inactive if value >= cutoff, else inconclusive
      "<" : MIC < value; active if value <= cutoff, else inconclusive
    """
    if operator == "=":
        return 1 if value <= cutoff else 0
    if operator == ">":
        return 0 if value >= cutoff else -1
    if operator == "<":
        return 1 if value <= cutoff else -1
    return -1


def aggregate_labels(labels):
    """Aggregate a list of (1, 0, -1) labels.

    Definitive results (0/1) take precedence; mean threshold at 0.5.
    Returns -1 if all inconclusive.
    """
    definitive = [l for l in labels if l != -1]
    if not definitive:
        return -1
    return 1 if (sum(definitive) / len(definitive)) >= 0.5 else 0


# ---------------------------------------------------------------------------
# Core binarisation: aggregate per std_smiles
# ---------------------------------------------------------------------------
def binarise(df):
    """Aggregate replicates per std_smiles and binarise at each cutoff.

    Returns DataFrame with columns:
        std_smiles, inchikey, mw, value, std, replicas, mic_10, mic_25, mic_50
    """
    if df.empty:
        return pd.DataFrame(columns=["std_smiles", "inchikey", "mw",
                                     "value", "std", "replicas"] + CUTOFF_COLS)

    rows = []
    for std_smiles, grp in df.groupby("std_smiles"):
        # Carry forward metadata from first row
        first = grp.iloc[0]
        inchikey = first["inchikey"]
        mw       = first["mw"]

        # Aggregate numeric values for reporting
        vals = grp["value"].dropna().tolist()
        n = len(vals)
        avg_val = round(float(np.mean(vals)), 3) if vals else np.nan
        std_val = round(float(np.std(vals, ddof=1)), 3) if n > 1 else 0.0

        # Aggregate operator: single value if uniform, "mixed" otherwise
        ops = grp["operator"].dropna().unique().tolist()
        agg_op = ops[0] if len(ops) == 1 else "mixed"

        # Binarise at each cutoff
        row = {
            "std_smiles": std_smiles,
            "inchikey":   inchikey,
            "mw":         round(float(mw), 3) if not pd.isna(mw) else np.nan,
            "value":      avg_val,
            "std":        std_val,
            "operator":   agg_op,
            "replicas":   n,
        }
        for cutoff, col in zip(CUTOFFS, CUTOFF_COLS):
            labels = [
                binarize_mic(r["operator"], r["value"], cutoff)
                for _, r in grp.iterrows()
                if not pd.isna(r["value"])
            ]
            row[col] = aggregate_labels(labels) if labels else -1

        rows.append(row)

    return pd.DataFrame(rows)[
        ["std_smiles", "inchikey", "mw", "value", "std", "operator", "replicas"] + CUTOFF_COLS
    ].reset_index(drop=True)


def save_dataset(df_bin, filepath, label):
    """Save dataset and return summary row."""
    n = len(df_bin)
    df_bin.to_csv(filepath, index=False)
    row = {"label": label, "n_compounds": n}
    for col in CUTOFF_COLS:
        definitive = df_bin[df_bin[col] != -1]
        n_active   = int((definitive[col] == 1).sum())
        n_definitive = len(definitive)
        row[f"n_active_{col}"]     = n_active
        row[f"n_inconclusive_{col}"] = int((df_bin[col] == -1).sum())
        row[f"active_rate_{col}"]  = round(n_active / n_definitive, 4) if n_definitive else np.nan
    print(f"  SAVED {label}: {n} compounds | "
          + " | ".join(f"{col}={row[f'active_rate_{col}']:.1%}"
                       for col in CUTOFF_COLS if not np.isnan(row[f"active_rate_{col}"])))
    return row


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
summary_rows = []

for patho_code in COADD_MIC_PATHO:
    fpath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    if not os.path.exists(fpath):
        print(f"\n[{patho_code}] Preprocessed file not found, skipping.")
        continue

    df = pd.read_csv(fpath)
    print(f"\n[{patho_code}] {len(df):,} rows loaded")

    patho_summary = []

    # Step 2: binarise per strain
    for strain_code in sorted(df["strain_code"].dropna().unique()):
        df_s   = df[df["strain_code"] == strain_code]
        df_bin = binarise(df_s)
        fname  = f"{patho_code}_{strain_code}.csv"
        row    = save_dataset(df_bin, os.path.join(binarised_dir, fname),
                              f"{patho_code}/{strain_code}")
        row["patho_code"]  = patho_code
        row["strain_code"] = strain_code
        row["file"]        = fname
        summary_rows.append(row)
        patho_summary.append({
            "label": strain_code,
            "n":     len(df_bin),
            **{col: (df_bin[df_bin[col] != -1][col].mean()
                     if (df_bin[col] != -1).any() else np.nan)
               for col in CUTOFF_COLS},
        })

    # Step 3: merged — only if more than one strain
    strains = df["strain_code"].dropna().unique()
    if len(strains) > 1:
        df_bin_merged = binarise(df)
        fname = f"{patho_code}_merged.csv"
        row   = save_dataset(df_bin_merged, os.path.join(binarised_dir, fname),
                             f"{patho_code}/merged")
        row["patho_code"]  = patho_code
        row["strain_code"] = "merged"
        row["file"]        = fname
        summary_rows.append(row)
        patho_summary.append({
            "label": "merged",
            "n":     len(df_bin_merged),
            **{col: (df_bin_merged[df_bin_merged[col] != -1][col].mean()
                     if (df_bin_merged[col] != -1).any() else np.nan)
               for col in CUTOFF_COLS},
        })
    else:
        print(f"  SKIP merged: only 1 strain present")

    # -----------------------------------------------------------------------
    # Step 4: plot
    # -----------------------------------------------------------------------
    if not patho_summary:
        continue

    import stylia

    # Format: print | Style: ersilia
    stylia.set_format("print")
    stylia.set_style("ersilia")

    labels   = [r["label"] for r in patho_summary]
    ns       = [r["n"] for r in patho_summary]
    rates    = {col: [r[col] if not pd.isna(r[col]) else 0
                      for r in patho_summary]
                for col in CUTOFF_COLS}
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
               label=f"≤{cutoff:.0f} µM")
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
    col_order = (
        ["patho_code", "strain_code", "n_compounds"]
        + [f"n_active_{col}"       for col in CUTOFF_COLS]
        + [f"n_inconclusive_{col}" for col in CUTOFF_COLS]
        + [f"active_rate_{col}"    for col in CUTOFF_COLS]
    )
    summary_df = pd.DataFrame(summary_rows)[col_order]
    summary_df.to_csv(os.path.join(output_dir, "summary.csv"), index=False)
    print(f"\nSummary written: {len(summary_df)} rows → output/05_binarise_coadd_mic/summary.csv")

print("\nDone.")
