"""
Summarise CO-ADD screening coverage from raw inhibition and dose-response files.

Reads config from config/manual/coadd_strains.csv.
Compound counts are unique SMILES only; no activity thresholds are applied.

Outputs (output/00_explore_coadd/):
  detailed_coverage_summary.csv  — one row per pathogen × strain_code × assay_type ×
                                   concentration; n_smiles = unique SMILES
  coverage_summary.csv           — one row per pathogen × strain_code × assay_type;
                                   n_smiles = unique SMILES across all concentrations
  compound_coverage.png          — bar chart: n_smiles per pathogen × assay type
  strain_overlap.png             — stacked horizontal bar per minority strain showing
                                   shared / majority-only / minority-only compound sets

Prints a SMILES comparison between the full CO-ADD library (inhibition + dose-response)
and the SPARK CO-ADD contribution file (data/raw/spark/SPARK Data CO-ADD Contribution.csv).
"""

import os
import sys
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

data_dir = os.path.join(root, "..", "data", "raw", "coadd")
outdir = os.path.join(root, "..", "output", "00_explore_coadd")
config_dir = os.path.join(root, "..", "config", "manual")

os.makedirs(outdir, exist_ok=True)

strains_cfg = pd.read_csv(os.path.join(config_dir, "coadd_strains.csv"))

# Build (organism, strain) → strain_code lookup
strain_lookup = {
    (row["ORGANISM"], row["STRAIN"]): row["strain_code"]
    for _, row in strains_cfg.iterrows()
}

# ---------------------------------------------------------------------------
# Inhibition data
# ---------------------------------------------------------------------------
print("Loading inhibition data...")
inhib = pd.read_csv(
    os.path.join(data_dir, "CO-ADD_InhibitionData_r03_01-02-2020_CSV.csv"),
    low_memory=False,
)
inhib = inhib[inhib["INHIB_AVE"].notnull()]
inhib = inhib[inhib["SMILES"].notnull()]

detailed_rows = []

for organism, org_strains in strains_cfg.groupby("ORGANISM"):
    code = org_strains.iloc[0]["patho_code"]
    df_p = inhib[inhib["ORGANISM"] == organism]
    if df_p.empty:
        continue

    print(f"  {code}: {len(df_p)} inhibition rows")
    print(f"    Strains: {df_p['STRAIN'].value_counts().to_dict()}")
    print(f"    Concentrations: {df_p['CONC'].value_counts().to_dict()}")

    # ALL strain: group by concentration
    for conc_val, df_c in df_p.groupby("CONC"):
        detailed_rows.append({
            "pathogen_code": code,
            "strain_code": "ALL",
            "assay_type": "inhibition",
            "concentration": conc_val,
            "n_smiles": df_c["SMILES"].nunique(),
        })

    # Per-strain breakdown
    for strain_val, df_s in df_p.groupby("STRAIN"):
        sc = strain_lookup.get((organism, strain_val), strain_val)
        for conc_val, df_sc in df_s.groupby("CONC"):
            detailed_rows.append({
                "pathogen_code": code,
                "strain_code": sc,
                "assay_type": "inhibition",
                "concentration": conc_val,
                "n_smiles": df_sc["SMILES"].nunique(),
            })

# ---------------------------------------------------------------------------
# Dose-response / MIC data
# ---------------------------------------------------------------------------
print("\nLoading dose-response data...")
dr = pd.read_csv(
    os.path.join(data_dir, "CO-ADD_DoseResponseData_r03_01-02-2020_CSV.csv"),
    low_memory=False,
)
dr = dr[dr["SMILES"].notnull()]
dr = dr[dr["ORGANISM"] != "Homo sapiens"]
dr_mic = dr[dr["DRVAL_TYPE"] == "MIC"].copy()
dr_mic = dr_mic[dr_mic["DRVAL_MEDIAN"].notnull()]
dr_mic = dr_mic[dr_mic["DRVAL_UNIT"].notnull()]

for organism, org_strains in strains_cfg.groupby("ORGANISM"):
    code = org_strains.iloc[0]["patho_code"]
    df_p = dr_mic[dr_mic["ORGANISM"] == organism]
    if df_p.empty:
        continue

    print(f"  {code}: {len(df_p)} MIC rows")
    print(f"    Strains: {df_p['STRAIN'].value_counts().to_dict()}")

    # ALL strain
    detailed_rows.append({
        "pathogen_code": code,
        "strain_code": "ALL",
        "assay_type": "mic",
        "concentration": "dose_response",
        "n_smiles": df_p["SMILES"].nunique(),
    })

    # Per-strain breakdown
    for strain_val, df_s in df_p.groupby("STRAIN"):
        sc = strain_lookup.get((organism, strain_val), strain_val)
        detailed_rows.append({
            "pathogen_code": code,
            "strain_code": sc,
            "assay_type": "mic",
            "concentration": "dose_response",
            "n_smiles": df_s["SMILES"].nunique(),
        })

# ---------------------------------------------------------------------------
# Save summaries
# ---------------------------------------------------------------------------
detailed = pd.DataFrame(detailed_rows).sort_values(
    ["pathogen_code", "strain_code", "assay_type", "concentration"]
).reset_index(drop=True)

detailed_path = os.path.join(outdir, "detailed_coverage_summary.csv")
detailed.to_csv(detailed_path, index=False)
print(f"\nDetailed coverage summary saved to {detailed_path}")

# Compressed: unique SMILES unioned across all concentrations per (pathogen, strain_code, assay_type)
# Built directly from raw dataframes so all tested SMILES are captured regardless of concentration group.
compressed_rows = []

for organism, org_strains in strains_cfg.groupby("ORGANISM"):
    code = org_strains.iloc[0]["patho_code"]

    # Inhibition
    df_p = inhib[inhib["ORGANISM"] == organism]
    if not df_p.empty:
        compressed_rows.append({
            "pathogen_code": code,
            "strain_code": "ALL",
            "assay_type": "inhibition",
            "n_smiles": df_p["SMILES"].nunique(),
        })
        for strain_val, df_s in df_p.groupby("STRAIN"):
            sc = strain_lookup.get((organism, strain_val), strain_val)
            compressed_rows.append({
                "pathogen_code": code,
                "strain_code": sc,
                "assay_type": "inhibition",
                "n_smiles": df_s["SMILES"].nunique(),
            })

    # MIC
    df_p = dr_mic[dr_mic["ORGANISM"] == organism]
    if not df_p.empty:
        compressed_rows.append({
            "pathogen_code": code,
            "strain_code": "ALL",
            "assay_type": "mic",
            "n_smiles": df_p["SMILES"].nunique(),
        })
        for strain_val, df_s in df_p.groupby("STRAIN"):
            sc = strain_lookup.get((organism, strain_val), strain_val)
            compressed_rows.append({
                "pathogen_code": code,
                "strain_code": sc,
                "assay_type": "mic",
                "n_smiles": df_s["SMILES"].nunique(),
            })

compressed = pd.DataFrame(compressed_rows).sort_values(
    ["pathogen_code", "strain_code", "assay_type"]
).reset_index(drop=True)

compressed_path = os.path.join(outdir, "coverage_summary.csv")
compressed.to_csv(compressed_path, index=False)
print(f"Coverage summary saved to {compressed_path}")

# Print strain coverage comparison
print("\n" + "=" * 60)
print("STRAIN COVERAGE: full library vs. subset")
print("=" * 60)

for (code, assay), grp in compressed.groupby(["pathogen_code", "assay_type"]):
    grp_strains = grp[grp["strain_code"] != "ALL"]
    if grp_strains.empty:
        continue
    max_n = grp_strains["n_smiles"].max()
    full = grp_strains[grp_strains["n_smiles"] == max_n]["strain_code"].tolist()
    subset = grp_strains[grp_strains["n_smiles"] < max_n][["strain_code", "n_smiles"]].sort_values(
        "n_smiles", ascending=False
    )
    print(f"\n  [{code}] {assay}  — full library ({int(max_n)} cpds): {full}")
    for _, r in subset.iterrows():
        pct = 100 * r["n_smiles"] / max_n
        print(f"    subset: {r['strain_code']}  {int(r['n_smiles'])} cpds  ({pct:.1f}%)")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
import stylia
stylia.set_format("print")
stylia.set_style("ersilia")

MIN_MOLECULES = 100

# Merged (all-strain) rows only
plot_df = compressed[compressed["strain_code"] == "ALL"].copy()
pathogens_order = plot_df["pathogen_code"].unique().tolist()
assay_types = ["inhibition", "mic"]

# ---------------------------------------------------------------------------
# Plot: Compound coverage bar chart
# n_smiles per pathogen × assay type (log scale, dashed MIN_MOLECULES)
# ---------------------------------------------------------------------------
print("\nPlot: Compound coverage...")

pal = stylia.CategoricalPalette("ersilia")
colors = pal.get(len(assay_types))
color_map = dict(zip(assay_types, colors))

n_pathogens = len(pathogens_order)
bar_width = 0.35
offsets = [-bar_width / 2, bar_width / 2]

fig, axs = stylia.create_figure(1, 1)
ax = axs.next()

for i, at in enumerate(assay_types):
    sub = plot_df[plot_df["assay_type"] == at].set_index("pathogen_code")
    xs = np.arange(n_pathogens)
    ys = [max(sub.loc[p, "n_smiles"], 1) if p in sub.index else 1 for p in pathogens_order]
    ax.bar(xs + offsets[i], ys, width=bar_width, color=color_map[at], label=at)

ax.axhline(MIN_MOLECULES, linestyle="--", color=stylia.NamedColors().gray, linewidth=1)
ax.set_yscale("log")
ax.set_xticks(np.arange(n_pathogens))
ax.set_xticklabels(pathogens_order, rotation=45, ha="right")
ax.legend(fontsize=stylia.FONTSIZE_SMALL)
stylia.label(ax, ylabel="Unique SMILES (log scale)", xlabel="",
             title="Compound coverage per pathogen  (dashed = 100 molecules)")
stylia.save_figure(os.path.join(outdir, "compound_coverage.png"))
print("  Saved compound_coverage.png")

# ---------------------------------------------------------------------------
# Plot: Compound overlap — majority vs minority strains (inhibition)
# For each organism with >1 strain, compare SMILES sets against the majority
# strain (most compounds). Shows: majority only | shared | minority only.
# ---------------------------------------------------------------------------
print("\nPlot: Strain overlap...")

overlap_rows = []
for organism, org_strains in strains_cfg.groupby("ORGANISM"):
    code = org_strains.iloc[0]["patho_code"]
    df_org = inhib[inhib["ORGANISM"] == organism]
    if df_org.empty:
        continue

    strain_smiles = {
        strain_lookup.get((organism, s), s): set(df_s["SMILES"].dropna())
        for s, df_s in df_org.groupby("STRAIN")
    }
    if len(strain_smiles) < 2:
        continue

    majority_strain = max(strain_smiles, key=lambda s: len(strain_smiles[s]))
    maj_set = strain_smiles[majority_strain]

    for strain, smi_set in strain_smiles.items():
        if strain == majority_strain:
            continue
        overlap_rows.append({
            "pathogen_code": code,
            "majority_strain": majority_strain,
            "minority_strain": strain,
            "shared": len(smi_set & maj_set),
            "majority_only": len(maj_set - smi_set),
            "minority_only": len(smi_set - maj_set),
        })

if overlap_rows:
    ov = pd.DataFrame(overlap_rows)
    labels = [f"{r['pathogen_code']}\n{r['minority_strain']}" for _, r in ov.iterrows()]
    n_bars = len(labels)

    pal = stylia.CategoricalPalette("ersilia")
    col_shared, col_maj_only, col_min_only = pal.get(3)

    fig, axs = stylia.create_figure(1, 1)
    ax = axs.next()

    ys = np.arange(n_bars)
    ax.barh(ys, ov["shared"],        color=col_shared,   label="shared with majority")
    ax.barh(ys, ov["majority_only"], color=col_maj_only, label="majority only",
            left=ov["shared"])
    ax.barh(ys, ov["minority_only"], color=col_min_only, label="minority only",
            left=ov["shared"] + ov["majority_only"])

    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.legend(fontsize=stylia.FONTSIZE_SMALL)
    stylia.label(ax, xlabel="Unique SMILES", ylabel="",
                 title="Compound overlap: minority strains vs majority strain (inhibition)")
    stylia.save_figure(os.path.join(outdir, "strain_overlap.png"))
    print("  Saved strain_overlap.png")
else:
    print("  No multi-strain organisms found, skipping.")

print("\nDone.")

# ---------------------------------------------------------------------------
# InChIKey comparison: full CO-ADD vs SPARK CO-ADD contribution
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("INCHIKEY COMPARISON: CO-ADD full library vs SPARK CO-ADD deposit")
print("=" * 60)

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi
from rdkit.Chem import rdinchi

def smiles_to_inchikey(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        inchi = MolToInchi(mol)
        if inchi is None:
            return None
        return rdinchi.InchiToInchiKey(inchi)
    except Exception:
        return None

print("  Computing InChIKeys for CO-ADD library...")
coadd_raw_smiles = set(inhib["SMILES"].dropna().unique()) | set(dr["SMILES"].dropna().unique())
coadd_keys = {k for smi in coadd_raw_smiles if (k := smiles_to_inchikey(smi)) is not None}

spark_coadd_path = os.path.join(root, "..", "data", "raw", "spark", "SPARK Data CO-ADD Contribution.csv")
spark_coadd = pd.read_csv(spark_coadd_path, low_memory=False)
print("  Computing InChIKeys for SPARK CO-ADD deposit...")
spark_raw_smiles = set(spark_coadd["SMILES"].dropna().unique())
spark_keys = {k for smi in spark_raw_smiles if (k := smiles_to_inchikey(smi)) is not None}

only_in_coadd = coadd_keys - spark_keys
only_in_spark = spark_keys - coadd_keys
shared = coadd_keys & spark_keys

print(f"\n  CO-ADD full library  : {len(coadd_keys):>8,} unique InChIKeys  (from {len(coadd_raw_smiles):,} SMILES)")
print(f"  SPARK CO-ADD deposit : {len(spark_keys):>8,} unique InChIKeys  (from {len(spark_raw_smiles):,} SMILES)")
print(f"  Shared               : {len(shared):>8,}")
print(f"  Only in CO-ADD       : {len(only_in_coadd):>8,}")
print(f"  Only in SPARK deposit: {len(only_in_spark):>8,}")
