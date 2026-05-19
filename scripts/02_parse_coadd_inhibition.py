"""
Parse CoADD inhibition (primary screening) data.

Step 1 – Preprocess raw data per pathogen:
  - Builds a lookup table of unique SMILES → std_smiles, inchikey, MW (computed once)
  - Extracts operator prefix from INHIB_AVE into a separate column
  - Converts assay concentration to µM using the dataset-wide average MW:
      * Already-µM entries (e.g. "25 uM") kept as-is
      * ug/mL entries: value * 1000 / avg_mw
      * Semicolon-separated concentrations (e.g. "32 ug/mL; 64 ug/mL") use
        the arithmetic mean of the stated values before conversion
  - Drops rows where std_smiles could not be computed or value is empty

Output: data/processed/coadd/02_preprocess_inh/{patho_code}.csv
Columns: std_smiles, inchikey, mw, patho_code, strain_code, conc, operator, value, std
"""

import os
import re
import sys
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.default import COADD_INH_PATHO
from src.utils import extract_operator

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
data_dir = os.path.join(root, "..", "data", "raw", "coadd")
preprocess_dir = os.path.join(root, "..", "data", "processed", "coadd", "02_preprocess_inh")
output_dir = os.path.join(root, "..", "output", "02_parse_coadd_inhibition")
config_dir = os.path.join(root, "..", "config", "manual")

os.makedirs(preprocess_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Load raw inhibition data
# ---------------------------------------------------------------------------
print("Loading CoADD inhibition data...")
df = pd.read_csv(
    os.path.join(data_dir, "CO-ADD_InhibitionData_r03_01-02-2020_CSV.csv"),
    low_memory=False,
)
df = df[df["SMILES"].notnull() & df["INHIB_AVE"].notnull()].copy()
print(f"Loaded {len(df):,} rows")

# ---------------------------------------------------------------------------
# Load strain config
# ---------------------------------------------------------------------------
strains_cfg = pd.read_csv(os.path.join(config_dir, "coadd_strains.csv"))
cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))
strain_map = dict(zip(strains_cfg["STRAIN"], strains_cfg["strain_code"]))
patho_organism = (
    strains_cfg[["patho_code", "ORGANISM"]]
    .drop_duplicates()
    .set_index("patho_code")["ORGANISM"]
    .to_dict()
)

# ---------------------------------------------------------------------------
# Load precomputed SMILES lookup (produced by 00_prepare_configs.py)
# ---------------------------------------------------------------------------
smiles_info_path = os.path.join(root, "..", "data", "processed", "coadd", "00_smiles_info.csv")
smiles_info = pd.read_csv(smiles_info_path)
smiles_lookup = smiles_info.set_index("smiles")[["std_smiles", "inchikey", "mw"]].to_dict("index")

avg_mw = float(smiles_info["mw"].dropna().mean())
print(f"Loaded {len(smiles_info):,} SMILES entries  |  Average MW: {avg_mw:.2f} g/mol")


# ---------------------------------------------------------------------------
# Helper: parse concentration string → µM
# ---------------------------------------------------------------------------
def parse_conc_to_um(conc_str, avg_mw):
    """Return assay concentration in µM.

    Handles:
      - "25 uM", "20 uM"           → value as-is
      - "32 ug/mL"                 → 32 * 1000 / avg_mw
      - "32 ug/mL; 64 ug/mL"      → mean(32, 64) * 1000 / avg_mw
    """
    if pd.isna(conc_str):
        return None
    s = str(conc_str).strip()
    nums = [float(x) for x in re.findall(r"[\d.]+", s)]
    if not nums:
        return None
    avg_num = float(np.mean(nums))
    if "uM" in s:
        return avg_num
    if "ug/mL" in s or "ug/ml" in s:
        return (avg_num / avg_mw) * 1000
    return None


# ---------------------------------------------------------------------------
# Step 1: Per-pathogen preprocessing
# ---------------------------------------------------------------------------
for patho_code in COADD_INH_PATHO:
    org_name = patho_organism.get(patho_code)
    if org_name is None:
        print(f"[{patho_code}] Not found in strain config, skipping.")
        continue

    df_p = df[df["ORGANISM"] == org_name].copy()
    if df_p.empty:
        print(f"[{patho_code}] No rows found, skipping.")
        continue

    print(f"\n[{patho_code}] {org_name} — {len(df_p):,} rows")

    # Map precomputed SMILES properties
    df_p["std_smiles"] = df_p["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("std_smiles"))
    df_p["inchikey"]   = df_p["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("inchikey"))
    df_p["mw"]         = df_p["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("mw"))

    # Operator + numeric value
    parsed = df_p["INHIB_AVE"].apply(lambda v: extract_operator(str(v)))
    df_p["operator"] = parsed.apply(lambda t: t[0])
    df_p["value"]    = parsed.apply(lambda t: float(t[1]) if t[1] is not None else None)

    # Concentration → µM
    df_p["concentration_um"] = df_p["CONC"].apply(lambda c: parse_conc_to_um(c, avg_mw))

    # Strain code and pathogen code
    df_p["strain_code"] = df_p["STRAIN"].map(strain_map)
    df_p["patho_code"]  = patho_code

    # Assemble, drop invalid rows, save
    out = (
        df_p[["std_smiles", "inchikey", "mw", "patho_code", "strain_code",
              "concentration_um", "operator", "value", "INHIB_STD"]]
        .rename(columns={"INHIB_STD": "std"})
        .dropna(subset=["std_smiles", "value"])
        .reset_index(drop=True)
    )

    filepath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    out.to_csv(filepath, index=False)
    print(f"  → {len(out):,} rows saved to {filepath}")

print("\nPreprocessing done.")

# ---------------------------------------------------------------------------
# Step 2: Distribution plots per pathogen
# ---------------------------------------------------------------------------
import stylia
from scipy.stats import gaussian_kde

# Format: print | Style: ersilia — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("ersilia")

inhib_cuts = cutoffs_cfg[cutoffs_cfg["assay_type"] == "inhib"]
CUTOFFS = sorted(set(
    inhib_cuts[["cutoff_low", "cutoff_mid", "cutoff_high"]].values.flatten().tolist()
))


def plot_conc_bar(ax, df):
    """Bar chart: measurement count per unique concentration (µM)."""
    counts = (
        df.dropna(subset=["concentration_um"])
        .groupby("concentration_um")["std_smiles"]
        .count()
        .reset_index()
        .sort_values("concentration_um")
    )
    nc = stylia.NamedColors()
    ax.bar(range(len(counts)), counts["std_smiles"], color=nc.blue)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels([f"{c:.1f}" for c in counts["concentration_um"]], rotation=45, ha="right")
    stylia.label(ax, xlabel="Concentration (µM)", ylabel="Measurements")


def plot_value_kde(ax, df):
    """KDE of inhibition % per strain; vertical dashed lines at cutoffs."""
    strains = sorted(df["strain_code"].dropna().unique())
    pal = stylia.CategoricalPalette("ersilia")
    colors = pal.get(max(len(strains), 1))
    vals_all = df["value"].dropna().values
    x_range = np.linspace(vals_all.min() - 5, vals_all.max() + 5, 500)
    for i, strain in enumerate(strains):
        vals = df[df["strain_code"] == strain]["value"].dropna().values
        if len(vals) < 10:
            continue
        ax.plot(x_range, gaussian_kde(vals)(x_range), color=colors[i], label=strain)
    nc = stylia.NamedColors()
    for cutoff in CUTOFFS:
        ax.axvline(cutoff, linestyle="--", color=nc.gray)
    ax.legend()
    stylia.label(ax, xlabel="Inhibition (%)", ylabel="Density")


print("\nGenerating distribution plots...")
for patho_code in COADD_INH_PATHO:
    fpath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    if not os.path.exists(fpath):
        print(f"  [{patho_code}] Preprocessed file not found, skipping.")
        continue
    df_p = pd.read_csv(fpath)

    fig, axs = stylia.create_figure(1, 2)
    plot_conc_bar(axs.next(), df_p)
    plot_value_kde(axs.next(), df_p)
    stylia.save_figure(os.path.join(output_dir, f"{patho_code}_distributions.png"))
    print(f"  → {patho_code}_distributions.png")

print("\nDone.")
