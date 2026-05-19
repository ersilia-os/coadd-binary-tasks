"""
Parse and binarise SPARK curated MIC data — all organisms, all strains.

Mimics COADD steps 04 (parse) + 05 (binarise), producing per-organism/strain files
in the same column format, using the curated MIC (µM) column from SPARK MIC Data.csv.

Outputs:
  data/processed/spark/10_preprocess_mic/{patho_code}.csv
    — one row per measurement: std_smiles, inchikey, mw, patho_code, strain_code, operator, value, std
  data/processed/spark/10_binarised_mic/{patho_code}_{strain_code}.csv
    — aggregated + binarised per strain: std_smiles, inchikey, mw, value, std, operator, replicas, mic_10, mic_25, mic_50
  data/processed/spark/10_binarised_mic/{patho_code}_merged.csv  (if >1 strain)
  output/10_parse_spark_mic/summary.csv
  output/10_parse_spark_mic/{patho_code}_mic_distribution.png
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchi
from rdkit.Chem import rdinchi
from scipy.stats import gaussian_kde

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.utils import extract_operator

try:
    from standardiser import standardise as std_lib
    HAS_STANDARDISER = True
except ImportError:
    HAS_STANDARDISER = False
    print("WARNING: standardiser not available — using RDKit canonical SMILES only")

# ─── paths ────────────────────────────────────────────────────────────────────

spark_mic_path = os.path.join(root, "..", "data", "raw", "spark", "SPARK MIC Data.csv")
preprocess_dir = os.path.join(root, "..", "data", "processed", "spark", "10_preprocess_mic")
binarised_dir  = os.path.join(root, "..", "data", "processed", "spark", "10_binarised_mic")
output_dir     = os.path.join(root, "..", "output", "10_parse_spark_mic")
config_dir     = os.path.join(root, "..", "config", "manual")

os.makedirs(preprocess_dir, exist_ok=True)
os.makedirs(binarised_dir,  exist_ok=True)
os.makedirs(output_dir,     exist_ok=True)

# ─── SPARK MIC column names ───────────────────────────────────────────────────

CUR_SPECIES = "Curated & Transformed MIC Data: Species"
CUR_MIC_UM  = "Curated & Transformed MIC Data: MIC (in µM) (µM)"
CUR_STRAIN  = "Curated & Transformed MIC Data: Strain"

# ─── cutoffs (reuse COADD config) ────────────────────────────────────────────

cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))
mic_cuts    = cutoffs_cfg[cutoffs_cfg["assay_type"] == "mic"]
CUTOFFS     = sorted(set(mic_cuts[["cutoff_low", "cutoff_mid", "cutoff_high"]].values.flatten()))
CUTOFF_COLS = [f"mic_{c:.0f}" for c in CUTOFFS]
print(f"Cutoffs: {CUTOFFS}  →  {CUTOFF_COLS}")

# ─── organism → patho_code mapping ───────────────────────────────────────────

strains_cfg = pd.read_csv(os.path.join(config_dir, "coadd_strains.csv"))
organism_to_patho = {
    row["ORGANISM"]: row["patho_code"]
    for _, row in strains_cfg.drop_duplicates("ORGANISM").iterrows()
}


def build_code_registry(species_list):
    """Build a collision-free species → patho_code mapping.

    Uses COADD config codes for known organisms; derives codes for the rest
    and appends a numeric suffix on collision.
    """
    registry  = {}
    used      = set(organism_to_patho.values())
    for name in sorted(species_list):  # alphabetical → deterministic
        if name in organism_to_patho:
            registry[name] = organism_to_patho[name]
            continue
        parts = name.split()
        base  = ((parts[0][0] + parts[1]).lower()
                 if len(parts) >= 2
                 else name.replace(" ", "").lower())
        code, i = base, 2
        while code in used:
            code = f"{base}{i}"
            i += 1
        registry[name] = code
        used.add(code)
    return registry


# ─── helpers ─────────────────────────────────────────────────────────────────

def parse_mic_um(val):
    """Extract (operator, float) from a SPARK curated µM value (string or number)."""
    if pd.isna(val):
        return None, None
    if isinstance(val, (int, float)):
        return "=", float(val)
    op, s = extract_operator(str(val))
    try:
        return op, float(s.strip())
    except (ValueError, AttributeError):
        return None, None


def normalise_strain(s):
    """'ATCC 25922' → 'ATCC25922'; NaN → 'unknown'."""
    if pd.isna(s):
        return "unknown"
    cleaned = re.sub(r"[^A-Za-z0-9]", "", str(s))
    return cleaned if cleaned else "unknown"


def standardise_smi(smi):
    """Return {std_smiles, inchikey, mw}; values are None on failure."""
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return {"std_smiles": None, "inchikey": None, "mw": None}
        if HAS_STANDARDISER:
            mol = std_lib.run(mol)
        std   = Chem.MolToSmiles(mol)
        inchi = MolToInchi(mol)
        key   = rdinchi.InchiToInchiKey(inchi) if inchi else None
        mw    = round(Descriptors.MolWt(mol), 3)
        return {"std_smiles": std, "inchikey": key, "mw": mw}
    except Exception:
        return {"std_smiles": None, "inchikey": None, "mw": None}


# ─── binarisation (identical logic to COADD step 05) ─────────────────────────

def binarize_mic(operator, value, cutoff):
    """1 = active, 0 = inactive, -1 = inconclusive. Lower MIC = more active."""
    if operator == "=":
        return 1 if value <= cutoff else 0
    if operator == ">":
        return 0 if value >= cutoff else -1
    if operator == "<":
        return 1 if value <= cutoff else -1
    return -1


def aggregate_labels(labels):
    definitive = [l for l in labels if l != -1]
    if not definitive:
        return -1
    return 1 if (sum(definitive) / len(definitive)) >= 0.5 else 0


def binarise(df):
    if df.empty:
        return pd.DataFrame(
            columns=["std_smiles", "inchikey", "mw",
                     "value", "std", "operator", "replicas"] + CUTOFF_COLS
        )
    rows = []
    for std_smiles, grp in df.groupby("std_smiles"):
        first  = grp.iloc[0]
        vals   = grp["value"].dropna().tolist()
        n      = len(vals)
        avg    = round(float(np.mean(vals)), 3) if vals else np.nan
        std_v  = round(float(np.std(vals, ddof=1)), 3) if n > 1 else 0.0
        ops    = grp["operator"].dropna().unique().tolist()
        agg_op = ops[0] if len(ops) == 1 else "mixed"

        row = {
            "std_smiles": std_smiles,
            "inchikey":   first["inchikey"],
            "mw":         round(float(first["mw"]), 3) if not pd.isna(first["mw"]) else np.nan,
            "value":      avg,
            "std":        std_v,
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


MIN_COMPOUNDS = 500  # minimum unique compounds to keep an individual-strain file


def save_dataset(df_bin, filepath, label, kept):
    """Compute summary stats and optionally save the file.

    filepath is written only when kept=True; otherwise only stats are returned.
    """
    n = len(df_bin)
    if kept:
        df_bin.to_csv(filepath, index=False)
    row = {"label": label, "n_compounds": n, "kept": kept}
    rate_parts = []
    for col in CUTOFF_COLS:
        definitive    = df_bin[df_bin[col] != -1]
        n_active      = int((definitive[col] == 1).sum())
        n_def         = len(definitive)
        active_rate   = round(n_active / n_def, 4) if n_def else np.nan
        row[f"n_active_{col}"]       = n_active
        row[f"n_inconclusive_{col}"] = int((df_bin[col] == -1).sum())
        row[f"active_rate_{col}"]    = active_rate
        if not (isinstance(active_rate, float) and np.isnan(active_rate)):
            rate_parts.append(f"{col}={active_rate:.1%}")
    tag = "KEPT" if kept else "skip"
    print(f"  [{tag}] {label}: {n} compounds | " + " | ".join(rate_parts))
    return row


# ─── load SPARK MIC data ──────────────────────────────────────────────────────

print("\nLoading SPARK MIC Data.csv...")
mic = pd.read_csv(spark_mic_path, low_memory=False)

curated = mic[
    mic[CUR_MIC_UM].notna() &
    mic["SMILES"].notna() &
    mic[CUR_SPECIES].notna()
].copy()
print(f"Curated MIC rows : {len(curated):,}")
print(f"Unique species   : {curated[CUR_SPECIES].nunique():,}")
print(f"Unique SMILES    : {curated['SMILES'].nunique():,}")

# ─── build SMILES lookup ──────────────────────────────────────────────────────

print("\nStandardising SMILES (this may take several minutes)...")
unique_smiles  = curated["SMILES"].dropna().unique()
smiles_lookup  = {}
for smi in tqdm(unique_smiles):
    smiles_lookup[smi] = standardise_smi(smi)

n_ok = sum(1 for v in smiles_lookup.values() if v["std_smiles"])
print(f"Lookup built: {n_ok:,} / {len(smiles_lookup):,} succeeded")

# ─── step 04: parse per organism ──────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 04: Parse per organism")
print("=" * 65)

species_list = sorted(curated[CUR_SPECIES].dropna().unique())
code_registry = build_code_registry(species_list)  # species → patho_code (collision-free)
patho_to_species = {}

for species in species_list:
    patho_code = code_registry[species]
    patho_to_species[patho_code] = species

    df_sp = curated[curated[CUR_SPECIES] == species].copy()

    df_sp["std_smiles"] = df_sp["SMILES"].map(lambda s: smiles_lookup.get(s, {}).get("std_smiles"))
    df_sp["inchikey"]   = df_sp["SMILES"].map(lambda s: smiles_lookup.get(s, {}).get("inchikey"))
    df_sp["mw"]         = df_sp["SMILES"].map(lambda s: smiles_lookup.get(s, {}).get("mw"))

    parsed           = df_sp[CUR_MIC_UM].apply(parse_mic_um)
    df_sp["operator"] = parsed.apply(lambda t: t[0])
    df_sp["value"]    = parsed.apply(lambda t: t[1])

    df_sp["strain_code"] = df_sp[CUR_STRAIN].apply(normalise_strain)
    df_sp["patho_code"]  = code_registry[species]
    df_sp["std"]         = np.nan

    out = (
        df_sp[["std_smiles", "inchikey", "mw", "patho_code", "strain_code",
               "operator", "value", "std"]]
        .dropna(subset=["std_smiles", "value"])
        .reset_index(drop=True)
    )

    out.to_csv(os.path.join(preprocess_dir, f"{patho_code}.csv"), index=False)
    ops     = out["operator"].value_counts().to_dict()
    n_strn  = out["strain_code"].nunique()
    print(f"  [{patho_code:<20}] {len(out):>6,} rows  |  {n_strn} strain(s)  |  ops: {ops}")


# ─── step 05: binarise per strain ─────────────────────────────────────────────

print("\n" + "=" * 65)
print("STEP 05: Binarise per strain")
print("=" * 65)

summary_rows = []

for patho_code in sorted(patho_to_species.keys()):
    fpath = os.path.join(preprocess_dir, f"{patho_code}.csv")
    if not os.path.exists(fpath):
        continue

    df = pd.read_csv(fpath)
    print(f"\n[{patho_code}] {len(df):,} rows")

    patho_summary = []

    for strain_code in sorted(df["strain_code"].dropna().unique()):
        df_s   = df[df["strain_code"] == strain_code]
        df_bin = binarise(df_s)
        fname  = f"{patho_code}_{strain_code}.csv"
        kept   = len(df_bin) >= MIN_COMPOUNDS
        row    = save_dataset(df_bin, os.path.join(binarised_dir, fname),
                              f"{patho_code}/{strain_code}", kept=kept)
        row.update({"patho_code": patho_code, "strain_code": strain_code, "file": fname})
        summary_rows.append(row)
        patho_summary.append({
            "label": strain_code,
            "n": len(df_bin),
            **{col: (df_bin[df_bin[col] != -1][col].mean()
                     if (df_bin[col] != -1).any() else np.nan)
               for col in CUTOFF_COLS},
        })

    strains = df["strain_code"].dropna().unique()
    if len(strains) > 1:
        df_merged = binarise(df)
        fname = f"{patho_code}_merged.csv"
        row   = save_dataset(df_merged, os.path.join(binarised_dir, fname),
                             f"{patho_code}/merged", kept=True)
        row.update({"patho_code": patho_code, "strain_code": "merged", "file": fname})
        summary_rows.append(row)
        patho_summary.append({
            "label": "merged",
            "n": len(df_merged),
            **{col: (df_merged[df_merged[col] != -1][col].mean()
                     if (df_merged[col] != -1).any() else np.nan)
               for col in CUTOFF_COLS},
        })

    # ── KDE plot ──────────────────────────────────────────────────────────────
    df_vals = df.dropna(subset=["value"])
    if len(df_vals) < 3:
        continue

    import stylia
    stylia.set_format("print")
    stylia.set_style("ersilia")

    nc = stylia.NamedColors()

    # Build strain list: ≥3 data points, capped at top-10 by compound count
    strain_counts = (
        df_vals.groupby("strain_code")["std_smiles"].nunique()
        .sort_values(ascending=False)
    )
    strains_with_data = [
        s for s in strain_counts.index
        if strain_counts[s] >= 3
    ][:10]

    if not strains_with_data:
        continue

    cm = stylia.CyclicColormap("ersilia")
    cm.fit(list(range(len(strains_with_data))))
    colors = cm.transform(list(range(len(strains_with_data))))

    all_vals = df_vals["value"].values
    x_min    = max(all_vals.min() * 0.5, 0.01)
    x_max    = all_vals.max() * 2
    x_range  = np.logspace(np.log10(x_min), np.log10(x_max), 500)

    fig, axs = stylia.create_figure(1, 1, width=0.5)
    ax = axs.next()

    for i, strain in enumerate(strains_with_data):
        vals = df_vals[df_vals["strain_code"] == strain]["value"].dropna().values
        kde  = gaussian_kde(np.log10(vals))
        ax.plot(x_range, kde(np.log10(x_range)), color=colors[i], label=strain)

    for cutoff in CUTOFFS:
        ax.axvline(cutoff, linestyle="--", color=nc.gray)

    ax.set_xscale("log")
    ax.legend()
    n_total_strains = len(strain_counts[strain_counts >= 3])
    suffix = f" (top 10 of {n_total_strains})" if n_total_strains > 10 else ""
    stylia.label(ax, xlabel="MIC (µM)", ylabel="Density", title=f"{patho_code}{suffix}")
    stylia.save_figure(os.path.join(output_dir, f"{patho_code}_mic_distribution.png"))
    print(f"  → {patho_code}_mic_distribution.png")


# ─── summary CSV ──────────────────────────────────────────────────────────────

if summary_rows:
    col_order = (
        ["patho_code", "strain_code", "n_compounds", "kept"]
        + [f"n_active_{col}"        for col in CUTOFF_COLS]
        + [f"n_inconclusive_{col}"  for col in CUTOFF_COLS]
        + [f"active_rate_{col}"     for col in CUTOFF_COLS]
    )
    summary_df = pd.DataFrame(summary_rows)[col_order]
    summary_df.to_csv(os.path.join(output_dir, "summary.csv"), index=False)
    print(f"\nSummary: {len(summary_df)} rows → output/10_parse_spark_mic/summary.csv")

print("\nDone.")
