"""
Parse CoADD MIC (dose-response) data.

Preprocessing per pathogen:
  - Loads MIC rows from the dose-response file
  - Maps precomputed std_smiles / inchikey / MW from 00_smiles_info.csv
  - Extracts censored operator from DRVAL_MEDIAN (e.g. ">10" → operator=">", value=10)
  - Converts MIC value to µM using molecule-specific MW (not average):
      * DRVAL_UNIT == "uM"    → kept as-is
      * DRVAL_UNIT == "ug/mL" → (numeric / mw) * 1000
  - Drops rows where std_smiles or value cannot be determined

Output: data/processed/coadd/04_preprocess_mic/{patho_code}.csv
Columns: std_smiles, inchikey, mw, patho_code, strain_code, operator, value, std

Note: there is no separate assay concentration for MIC — value IS the MIC in µM.
      std is NaN (no replicate std available in the raw file).
"""

import os
import sys
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.default import COADD_MIC_PATHO
from src.utils import extract_operator

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
data_dir       = os.path.join(root, "..", "data", "raw", "coadd")
preprocess_dir = os.path.join(root, "..", "data", "processed", "coadd", "04_preprocess_mic")
output_dir     = os.path.join(root, "..", "output", "04_parse_coadd_mic")
config_dir     = os.path.join(root, "..", "config", "manual")

os.makedirs(preprocess_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Load raw MIC data
# ---------------------------------------------------------------------------
print("Loading CoADD dose-response data...")
df = pd.read_csv(
    os.path.join(data_dir, "CO-ADD_DoseResponseData_r03_01-02-2020_CSV.csv"),
    low_memory=False,
)
df = df[
    (df["DRVAL_TYPE"] == "MIC") &
    df["SMILES"].notnull() &
    df["DRVAL_MEDIAN"].notnull()
].copy()
print(f"Loaded {len(df):,} MIC rows")
print(f"Units: {df['DRVAL_UNIT'].value_counts().to_dict()}")

# ---------------------------------------------------------------------------
# Load strain config and SMILES lookup
# ---------------------------------------------------------------------------
strains_cfg = pd.read_csv(os.path.join(config_dir, "coadd_strains.csv"))
strain_map = dict(zip(strains_cfg["STRAIN"], strains_cfg["strain_code"]))
patho_organism = (
    strains_cfg[["patho_code", "ORGANISM"]]
    .drop_duplicates()
    .set_index("patho_code")["ORGANISM"]
    .to_dict()
)

smiles_info = pd.read_csv(
    os.path.join(root, "..", "data", "processed", "coadd", "00_smiles_info.csv")
)
smiles_lookup = smiles_info.set_index("smiles")[["std_smiles", "inchikey", "mw"]].to_dict("index")
print(f"Loaded {len(smiles_info):,} SMILES entries from lookup")


# ---------------------------------------------------------------------------
# Helper: convert MIC value to µM using molecule-specific MW
# ---------------------------------------------------------------------------
def mic_to_um(numeric, unit, mw):
    """Return MIC in µM. Uses molecule-specific MW for ug/mL conversion."""
    if pd.isna(numeric) or pd.isna(mw) or mw <= 0:
        return None
    if unit == "uM":
        return float(numeric)
    if unit == "ug/mL":
        return (float(numeric) / mw) * 1000
    return None


# ---------------------------------------------------------------------------
# Per-pathogen preprocessing
# ---------------------------------------------------------------------------
for patho_code in COADD_MIC_PATHO:
    org_name = patho_organism.get(patho_code)
    if org_name is None:
        print(f"\n[{patho_code}] Not found in strain config, skipping.")
        continue

    df_p = df[df["ORGANISM"] == org_name].copy()
    if df_p.empty:
        print(f"\n[{patho_code}] No MIC rows found, skipping.")
        continue

    print(f"\n[{patho_code}] {org_name} — {len(df_p):,} rows")

    # Map precomputed SMILES properties
    df_p["std_smiles"] = df_p["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("std_smiles"))
    df_p["inchikey"]   = df_p["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("inchikey"))
    df_p["mw"]         = df_p["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("mw"))
    df_p["mw"]         = df_p["mw"].round(3)

    # Extract operator + numeric from DRVAL_MEDIAN
    parsed = df_p["DRVAL_MEDIAN"].apply(lambda v: extract_operator(str(v)))
    df_p["operator"] = parsed.apply(lambda t: t[0])
    df_p["numeric"]  = parsed.apply(lambda t: float(t[1]) if t[1] is not None else None)

    # Convert to µM using molecule-specific MW
    df_p["value"] = df_p.apply(
        lambda r: mic_to_um(r["numeric"], r["DRVAL_UNIT"], r["mw"]), axis=1
    )

    # Strain and pathogen codes; std not available in raw MIC data
    df_p["strain_code"] = df_p["STRAIN"].map(strain_map)
    df_p["patho_code"]  = patho_code
    df_p["std"]         = np.nan

    # Assemble, drop invalid rows, save
    out = (
        df_p[["std_smiles", "inchikey", "mw", "patho_code", "strain_code",
              "operator", "value", "std"]]
        .dropna(subset=["std_smiles", "value"])
        .reset_index(drop=True)
    )

    filepath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    out.to_csv(filepath, index=False)
    op_counts = out["operator"].value_counts().to_dict()
    print(f"  → {len(out):,} rows saved | operators: {op_counts}")

print("\nPreprocessing done.")

# ---------------------------------------------------------------------------
# KDE plots: MIC value distribution per pathogen, one line per strain
# ---------------------------------------------------------------------------
import stylia
from scipy.stats import gaussian_kde

# Format: print | Style: ersilia
stylia.set_format("print")
stylia.set_style("ersilia")

cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))
mic_cuts = cutoffs_cfg[cutoffs_cfg["assay_type"] == "mic"]
CUTOFFS = sorted(set(
    mic_cuts[["cutoff_low", "cutoff_mid", "cutoff_high"]].values.flatten().tolist()
))

print("\nGenerating MIC distribution plots...")
for patho_code in COADD_MIC_PATHO:
    fpath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    if not os.path.exists(fpath):
        continue
    df_p = pd.read_csv(fpath)

    df_vals = df_p.dropna(subset=["value"])
    if df_vals.empty:
        print(f"  [{patho_code}] No values for KDE, skipping.")
        continue

    strains = sorted(df_vals["strain_code"].dropna().unique())
    pal    = stylia.CategoricalPalette("ersilia")
    colors = pal.get(max(len(strains), 1))
    nc     = stylia.NamedColors()

    all_vals = df_vals["value"].values
    x_min = max(all_vals.min() * 0.5, 0.01)
    x_max = all_vals.max() * 2
    x_range = np.logspace(np.log10(x_min), np.log10(x_max), 500)

    fig, axs = stylia.create_figure(1, 1, width=0.5)
    ax = axs.next()

    for i, strain in enumerate(strains):
        vals = df_vals[df_vals["strain_code"] == strain]["value"].dropna().values
        if len(vals) < 3:
            continue
        log_vals = np.log10(vals)
        kde = gaussian_kde(log_vals)
        ax.plot(x_range, kde(np.log10(x_range)), color=colors[i], label=strain)

    for cutoff in CUTOFFS:
        ax.axvline(cutoff, linestyle="--", color=nc.gray)

    ax.set_xscale("log")
    ax.legend()
    stylia.label(ax, xlabel="MIC (µM)", ylabel="Density", title=patho_code)

    stylia.save_figure(os.path.join(output_dir, f"{patho_code}_mic_distribution.png"))
    print(f"  → {patho_code}_mic_distribution.png")

print("\nDone.")
