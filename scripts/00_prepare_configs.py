"""
Generate config skeleton files from raw CO-ADD and SPARK data.

Extracts available organisms/strains and assay types directly from raw data files,
with compound counts. The resulting files are intended for manual review and editing
(adding pathogen codes, strain codes, binarization cutoffs, etc.).

Output files:
  config/coadd_strains.csv   — unique ORGANISM × STRAIN from CO-ADD raw files
  config/coadd_cutoffs.csv   — concentration groups and value units from CO-ADD
  config/spark_strains.csv   — unique species × strain per SPARK file and assay type
  config/spark_cutoffs.csv   — available value columns per SPARK file and assay type
"""

import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from standardiser import standardise
from tqdm import tqdm

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

coadd_dir = os.path.join(root, "..", "data", "raw", "coadd")
spark_dir = os.path.join(root, "..", "data", "raw", "spark")
config_dir = os.path.join(root, "..", "config")
processed_dir = os.path.join(root, "..", "data", "processed", "coadd")

os.makedirs(processed_dir, exist_ok=True)


# =============================================================================
# CO-ADD
# =============================================================================
print("=" * 60)
print("CO-ADD")
print("=" * 60)

# --- Inhibition (primary screening) ---
print("\nLoading CO-ADD inhibition data...")
inhib = pd.read_csv(
    os.path.join(coadd_dir, "CO-ADD_InhibitionData_r03_01-02-2020_CSV.csv"),
    usecols=["ORGANISM", "STRAIN", "SMILES", "CONC"],
    low_memory=False,
)
inhib = inhib[inhib["SMILES"].notnull()]
print(f"  {len(inhib)} rows, {inhib['SMILES'].nunique()} unique SMILES")

# --- Dose-response (secondary screening, MIC only) ---
print("Loading CO-ADD dose-response data...")
dr = pd.read_csv(
    os.path.join(coadd_dir, "CO-ADD_DoseResponseData_r03_01-02-2020_CSV.csv"),
    usecols=["ORGANISM", "STRAIN", "SMILES", "DRVAL_TYPE", "DRVAL_UNIT"],
    low_memory=False,
)
dr = dr[dr["SMILES"].notnull()]
dr = dr[dr["ORGANISM"] != "Homo sapiens"]
dr_mic = dr[dr["DRVAL_TYPE"] == "MIC"]
print(f"  {len(dr_mic)} MIC rows, {dr_mic['SMILES'].nunique()} unique SMILES")

# ---------------------------------------------------------------------------
# 00_smiles_info.csv
# Unique SMILES (from both inhibition and MIC files) → std_smiles, inchikey, mw
# ---------------------------------------------------------------------------
smiles_info_path = os.path.join(processed_dir, "00_smiles_info.csv")

if os.path.exists(smiles_info_path):
    print(f"\n00_smiles_info.csv already exists, skipping computation.")
else:
    all_smiles = pd.Series(
        pd.concat([inhib["SMILES"], dr["SMILES"]]).dropna().unique(),
        name="smiles",
    )
    print(f"\nProcessing {len(all_smiles):,} unique SMILES (inhib + MIC)...")

    rows = []
    for smi in tqdm(all_smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                rows.append({"smiles": smi, "std_smiles": None, "inchikey": None, "mw": None})
                continue
            mol = standardise.run(mol)
            rows.append({
                "smiles": smi,
                "std_smiles": Chem.MolToSmiles(mol),
                "inchikey": InchiToInchiKey(MolToInchi(mol)),
                "mw": Descriptors.MolWt(mol),
            })
        except Exception:
            rows.append({"smiles": smi, "std_smiles": None, "inchikey": None, "mw": None})

    smiles_info = pd.DataFrame(rows)
    smiles_info.to_csv(smiles_info_path, index=False)
    n_valid = smiles_info["std_smiles"].notnull().sum()
    print(f"00_smiles_info.csv written: {n_valid:,} valid / {len(smiles_info):,} total")

# ---------------------------------------------------------------------------
# coadd_strains.csv
# One row per ORGANISM × STRAIN, with compound counts for each assay type.
# ---------------------------------------------------------------------------
inhib_grp = (
    inhib.groupby(["ORGANISM", "STRAIN"])["SMILES"]
    .nunique()
    .reset_index()
    .rename(columns={"SMILES": "inhib_n_compounds"})
)
mic_grp = (
    dr_mic.groupby(["ORGANISM", "STRAIN"])["SMILES"]
    .nunique()
    .reset_index()
    .rename(columns={"SMILES": "mic_n_compounds"})
)
coadd_strains = (
    inhib_grp.merge(mic_grp, on=["ORGANISM", "STRAIN"], how="outer")
    .sort_values(["ORGANISM", "STRAIN"])
    .reset_index(drop=True)
)
coadd_strains.to_csv(os.path.join(config_dir, "coadd_strains.csv"), index=False)
print(f"\ncoadd_strains.csv written: {len(coadd_strains)} rows")

# ---------------------------------------------------------------------------
# coadd_cutoffs.csv
# For inhibition: one row per CONC value with compound count.
# For MIC: one row per unit type with compound count.
# ---------------------------------------------------------------------------
inhib_concs = (
    inhib.groupby("CONC")["SMILES"]
    .nunique()
    .reset_index()
    .rename(columns={"CONC": "value_group", "SMILES": "n_compounds"})
)
inhib_concs["assay_type"] = "inhib"
inhib_concs["unit"] = "%"

mic_units = (
    dr_mic.groupby("DRVAL_UNIT")["SMILES"]
    .nunique()
    .reset_index()
    .rename(columns={"DRVAL_UNIT": "value_group", "SMILES": "n_compounds"})
)
mic_units["assay_type"] = "mic"
mic_units["unit"] = mic_units["value_group"]

coadd_cutoffs = pd.concat(
    [
        inhib_concs[["assay_type", "value_group", "unit", "n_compounds"]],
        mic_units[["assay_type", "value_group", "unit", "n_compounds"]],
    ],
    ignore_index=True,
)
coadd_cutoffs.to_csv(os.path.join(config_dir, "coadd_cutoffs.csv"), index=False)
print(f"coadd_cutoffs.csv written: {len(coadd_cutoffs)} rows")


# =============================================================================
# SPARK
# =============================================================================
print("\n" + "=" * 60)
print("SPARK")
print("=" * 60)

# Files that contain activity + species/strain data.
# "SPARK Data Compounds & Physicochemical Properties.csv" is excluded (no activity).
SPARK_FILES = [
    ("SPARK MIC Data.csv", "mic"),
    ("SPARK IC50 Data.csv", "ic50"),
    ("SPARK Accumulation Data.csv", "accumulation"),
    ("SPARK Data CO-ADD Contribution.csv", "co-add"),
    ("SPARK Data Achaogen Contribution.csv", "achaogen"),
    ("SPARK Data Merck & Kyorin Contribution.csv", "merck_kyorin"),
    ("SPARK Data Novartis Contribution.csv", "novartis"),
    ("SPARK Data Quave Lab {Emory University} Publications.csv", "quave"),
]

# Keywords that identify value columns worth reporting in spark_cutoffs.csv
VALUE_KEYWORDS = ["mic", "ic50", "inhibition %", "inhibition perc", "accumulated compound"]


def find_species_strain_pairs(columns):
    """
    Return list of (species_col, strain_col) pairs that share the same prefix
    (i.e. the text before the last ': Species' / ': Strain').
    """
    pairs = []
    for col in columns:
        if col.endswith("Species") or col.endswith(" Species"):
            # Try to find a matching Strain column with the same prefix
            prefix = col[: col.rfind("Species")]
            strain_col = prefix + "Strain"
            if strain_col in columns:
                pairs.append((col, strain_col))
    # Also handle bare "Species" / "Strain" columns (e.g. Accumulation file)
    if "Species" in columns and "Strain" in columns:
        if ("Species", "Strain") not in pairs:
            pairs.append(("Species", "Strain"))
    return pairs


def infer_assay_type(col_name):
    col_lower = col_name.lower()
    if "mic" in col_lower:
        return "mic"
    if "ic50" in col_lower:
        return "ic50"
    if "inhibit" in col_lower:
        return "inhibition"
    if "accum" in col_lower:
        return "accumulation"
    return "unknown"


def infer_unit(col_name):
    for unit in ["µM", "uM", "nM", "µg/mL", "ug/mL", "nmol", "pg", "ng/mL", "%"]:
        if unit in col_name:
            return unit
    return ""


spark_strains_rows = []
spark_cutoffs_rows = []

for fname, source_key in SPARK_FILES:
    fpath = os.path.join(spark_dir, fname)
    if not os.path.exists(fpath):
        print(f"\n  SKIP {fname}: file not found")
        continue

    print(f"\nProcessing: {fname}")

    # Read column names only first to avoid loading huge files unnecessarily
    header_df = pd.read_csv(fpath, nrows=0, low_memory=False)
    all_cols = list(header_df.columns)

    # Find SMILES column
    smiles_col = next((c for c in all_cols if c.strip().lower() == "smiles"), None)

    # Find (species, strain) column pairs
    pairs = find_species_strain_pairs(all_cols)
    print(f"  SMILES col: {smiles_col}")
    print(f"  Species/Strain pairs: {pairs}")

    # Find relevant value columns
    value_cols = [
        c for c in all_cols
        if any(kw in c.lower() for kw in VALUE_KEYWORDS)
        and "extracted" not in c.lower()       # prefer curated over extracted
        and "terms of use" not in c.lower()
    ]
    print(f"  Value cols: {value_cols}")

    if not pairs:
        print("  No species/strain columns found, skipping.")
        continue

    # Determine which columns to actually load
    load_cols = set()
    if smiles_col:
        load_cols.add(smiles_col)
    for sp_col, st_col in pairs:
        load_cols.add(sp_col)
        load_cols.add(st_col)

    df = pd.read_csv(fpath, usecols=list(load_cols), low_memory=False)
    if smiles_col:
        df = df[df[smiles_col].notnull()]

    # Extract unique species × strain per pair
    for sp_col, st_col in pairs:
        # Infer assay type from the column name prefix
        assay_type = infer_assay_type(sp_col)

        if smiles_col:
            grp = (
                df.dropna(subset=[sp_col, st_col])
                .groupby([sp_col, st_col])[smiles_col]
                .nunique()
                .reset_index()
                .rename(columns={sp_col: "species", st_col: "strain", smiles_col: "n_compounds"})
            )
        else:
            grp = (
                df.dropna(subset=[sp_col, st_col])
                .groupby([sp_col, st_col])
                .size()
                .reset_index(name="n_compounds")
                .rename(columns={sp_col: "species", st_col: "strain"})
            )

        grp["source_file"] = source_key
        grp["assay_type"] = assay_type
        spark_strains_rows.append(grp)
        print(f"    {sp_col}: {len(grp)} unique (species, strain) pairs")

    # Record available value columns (no data loading needed for this)
    for val_col in value_cols:
        spark_cutoffs_rows.append({
            "source_file": source_key,
            "assay_type": infer_assay_type(val_col),
            "value_column": val_col,
            "unit": infer_unit(val_col),
        })

# ---------------------------------------------------------------------------
# spark_strains.csv
# ---------------------------------------------------------------------------
if spark_strains_rows:
    spark_strains = (
        pd.concat(spark_strains_rows, ignore_index=True)
        .sort_values(["source_file", "assay_type", "species", "strain"])
        .reset_index(drop=True)
    )
    spark_strains = spark_strains[["source_file", "assay_type", "species", "strain", "n_compounds"]]
    spark_strains.to_csv(os.path.join(config_dir, "spark_strains.csv"), index=False)
    print(f"\nspark_strains.csv written: {len(spark_strains)} rows")
else:
    print("\nNo SPARK strain data extracted.")

# ---------------------------------------------------------------------------
# spark_cutoffs.csv
# ---------------------------------------------------------------------------
if spark_cutoffs_rows:
    spark_cutoffs = pd.DataFrame(spark_cutoffs_rows)[
        ["source_file", "assay_type", "value_column", "unit"]
    ]
    spark_cutoffs.to_csv(os.path.join(config_dir, "spark_cutoffs.csv"), index=False)
    print(f"spark_cutoffs.csv written: {len(spark_cutoffs)} rows")
else:
    print("No SPARK cutoffs data extracted.")

print("\nDone.")
