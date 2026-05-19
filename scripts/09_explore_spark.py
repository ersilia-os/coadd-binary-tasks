"""
Independent overview of SPARK data.

Answers:
  1. Are lab-contribution files (Achaogen, CO-ADD, Merck & Kyorin, Novartis, Quave)
     already captured in the three summary assay files (MIC, IC50, Accumulation)?
  2. Which organisms are tested, in which assays, and with how many unique compounds?
     (No pipeline filter — all species reported.)
  3. What is the compound overlap between SPARK and the already-processed CO-ADD dataset?

Outputs (output/09_explore_spark/):
  file_summary.csv             — one row per SPARK file
  contribution_vs_summary.csv  — coverage of each lab's SMILES in summary files
  organism_coverage.csv        — unique SMILES per species × assay type
  spark_vs_coadd.csv           — per-species InChIKey overlap with CO-ADD processed
  spark_mic_coverage_all.png   — top-25 species by curated MIC compound count
  spark_vs_coadd_overlap.png   — top-20 species SPARK vs CO-ADD overlap
"""

import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi
from rdkit.Chem import rdinchi

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.default import COADD_MIC_PATHO, COADD_INH_PATHO

spark_dir    = os.path.join(root, "..", "data", "raw", "spark")
processed_dir = os.path.join(root, "..", "data", "processed", "coadd")
config_dir   = os.path.join(root, "..", "config", "manual")
outdir       = os.path.join(root, "..", "output", "09_explore_spark")
os.makedirs(outdir, exist_ok=True)


# ─── helpers ─────────────────────────────────────────────────────────────────

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


def inchi_to_inchikey(inchi):
    try:
        if not isinstance(inchi, str):
            return None
        key = rdinchi.InchiToInchiKey(inchi)
        return key if key else None
    except Exception:
        return None


def smiles_set_to_inchikeys(smiles_iter):
    keys = set()
    for smi in smiles_iter:
        k = smiles_to_inchikey(smi)
        if k:
            keys.add(k)
    return keys


# ─── config ──────────────────────────────────────────────────────────────────

strains_cfg = pd.read_csv(os.path.join(config_dir, "coadd_strains.csv"))
organism_to_patho = {
    row["ORGANISM"]: row["patho_code"]
    for _, row in strains_cfg.drop_duplicates("ORGANISM").iterrows()
}
pipeline_pathos = set(COADD_MIC_PATHO) | set(COADD_INH_PATHO)


# ─── file catalogue ──────────────────────────────────────────────────────────

SPARK_FILES = [
    ("SPARK Accumulation Data.csv",                              "accumulation"),
    ("SPARK Data Achaogen Contribution.csv",                     "mic+ic50+pk"),
    ("SPARK Data CO-ADD Contribution.csv",                       "mic+inhibition"),
    ("SPARK Data Compounds & Physicochemical Properties.csv",    "registry+pchem"),
    ("SPARK Data Merck & Kyorin Contribution.csv",               "mic"),
    ("SPARK Data Novartis Contribution.csv",                     "mic+ic50"),
    ("SPARK Data Quave Lab {Emory University} Publications.csv", "mic"),
    ("SPARK IC50 Data.csv",                                      "ic50"),
    ("SPARK MIC Data.csv",                                       "mic"),
]

SUMMARY_FILES = {
    "SPARK MIC Data.csv",
    "SPARK IC50 Data.csv",
    "SPARK Accumulation Data.csv",
}

CONTRIBUTOR_FILES = {
    "SPARK Data Achaogen Contribution.csv":                     "Achaogen",
    "SPARK Data CO-ADD Contribution.csv":                       "CO-ADD",
    "SPARK Data Merck & Kyorin Contribution.csv":               "Merck & Kyorin",
    "SPARK Data Novartis Contribution.csv":                     "Novartis",
    "SPARK Data Quave Lab {Emory University} Publications.csv": "Quave Lab",
}


# ─── section 1: per-file stats ────────────────────────────────────────────────

print("=" * 65)
print("SECTION 1: Per-file statistics")
print("=" * 65)

cached_dfs        = {}  # fname -> DataFrame  (summary files only)
cached_contrib_smiles = {}  # fname -> set of SMILES  (contributor files)

file_rows = []
for fname, assay_type in SPARK_FILES:
    fpath = os.path.join(spark_dir, fname)
    df = pd.read_csv(fpath, low_memory=False)

    if fname in SUMMARY_FILES:
        cached_dfs[fname] = df
    if fname in CONTRIBUTOR_FILES:
        cached_contrib_smiles[fname] = set(df["SMILES"].dropna().unique()) if "SMILES" in df.columns else set()

    n_rows   = len(df)
    n_smiles = df["SMILES"].dropna().nunique() if "SMILES" in df.columns else 0

    if "InChI" in df.columns:
        keys = {inchi_to_inchikey(i) for i in df["InChI"].dropna().unique()}
        keys.discard(None)
        n_inchikeys = len(keys)
        key_note = "(from InChI col)"
    else:
        unique_smiles = df["SMILES"].dropna().unique() if "SMILES" in df.columns else []
        keys = smiles_set_to_inchikeys(unique_smiles)
        n_inchikeys = len(keys)
        key_note = "(from SMILES)"

    file_rows.append({
        "file": fname,
        "assay_type": assay_type,
        "n_rows": n_rows,
        "n_smiles": n_smiles,
        "n_inchikeys": n_inchikeys,
    })
    print(f"  {fname}")
    print(f"    {n_rows:>10,} rows  |  {n_smiles:>8,} SMILES  |  {n_inchikeys:>8,} InChIKeys {key_note}")

file_summary = pd.DataFrame(file_rows)
file_summary.to_csv(os.path.join(outdir, "file_summary.csv"), index=False)
print(f"\nSaved file_summary.csv")


# ─── section 2: are lab contributions in the summary files? ──────────────────

print("\n" + "=" * 65)
print("SECTION 2: Lab contribution files vs summary files")
print("=" * 65)
print("""
  Checks what fraction of each lab's SMILES appear in the three
  SPARK summary files (MIC, IC50, Accumulation).

  Note: CO-ADD's single-concentration inhibition data is NOT stored
  in any summary file; only its MIC records overlap with SPARK MIC Data.
""")

mic_smiles = set(cached_dfs["SPARK MIC Data.csv"]["SMILES"].dropna().unique())
ic50_smiles = set(cached_dfs["SPARK IC50 Data.csv"]["SMILES"].dropna().unique())
acc_smiles  = set(cached_dfs["SPARK Accumulation Data.csv"]["SMILES"].dropna().unique())

print(f"  Summary file SMILES counts:")
print(f"    SPARK MIC Data         : {len(mic_smiles):>7,}")
print(f"    SPARK IC50 Data        : {len(ic50_smiles):>7,}")
print(f"    SPARK Accumulation Data: {len(acc_smiles):>7,}")
print()

contrib_rows = []
for fname, label in CONTRIBUTOR_FILES.items():
    smiles = cached_contrib_smiles[fname]
    n = len(smiles)
    in_mic  = smiles & mic_smiles
    in_ic50 = smiles & ic50_smiles
    in_acc  = smiles & acc_smiles
    in_any  = in_mic | in_ic50 | in_acc
    not_in_any = smiles - in_any
    contrib_rows.append({
        "contributor": label,
        "n_smiles": n,
        "n_in_mic": len(in_mic),
        "pct_in_mic": round(100 * len(in_mic) / n, 1) if n else 0.0,
        "n_in_ic50": len(in_ic50),
        "pct_in_ic50": round(100 * len(in_ic50) / n, 1) if n else 0.0,
        "n_in_accumulation": len(in_acc),
        "pct_in_accumulation": round(100 * len(in_acc) / n, 1) if n else 0.0,
        "n_not_in_any": len(not_in_any),
        "pct_not_in_any": round(100 * len(not_in_any) / n, 1) if n else 0.0,
    })
    pct_mic  = 100 * len(in_mic)  / n if n else 0
    pct_ic50 = 100 * len(in_ic50) / n if n else 0
    pct_acc  = 100 * len(in_acc)  / n if n else 0
    pct_none = 100 * len(not_in_any) / n if n else 0
    print(f"  {label:<20}  {n:>6,} SMILES")
    print(f"    → in MIC       : {len(in_mic):>6,}  ({pct_mic:.1f}%)")
    print(f"    → in IC50      : {len(in_ic50):>6,}  ({pct_ic50:.1f}%)")
    print(f"    → in Accum     : {len(in_acc):>6,}  ({pct_acc:.1f}%)")
    print(f"    → not in any   : {len(not_in_any):>6,}  ({pct_none:.1f}%)")

contrib_df = pd.DataFrame(contrib_rows)
contrib_df.to_csv(os.path.join(outdir, "contribution_vs_summary.csv"), index=False)
print(f"\nSaved contribution_vs_summary.csv")


# ─── section 3: comprehensive organism × assay coverage ──────────────────────

print("\n" + "=" * 65)
print("SECTION 3: Organism × assay coverage (all species, no pipeline filter)")
print("=" * 65)

CUR_SPECIES  = "Curated & Transformed MIC Data: Species"
CUR_MIC_VAL  = "Curated & Transformed MIC Data: MIC (in µM) (µM)"
EXT_SPECIES  = "Extracted & Uploaded MIC Data: Species"
EXT_MIC_VAL  = "Extracted & Uploaded MIC Data: MIC (in µM) (µM)"
IC50_SPECIES = "Curated & Transformed IC50 Data: Species"
IC50_VAL     = "Curated & Transformed IC50 Data: IC50 (uM) (uM)"
ACC_SPECIES  = "Curated & Transformed Accumulation Data: Species"
ACC_VAL      = "Curated & Transformed Accumulation Data: Accumulated Compound (μM)"

mic    = cached_dfs["SPARK MIC Data.csv"]
ic50_df = cached_dfs["SPARK IC50 Data.csv"]
acc_df  = cached_dfs["SPARK Accumulation Data.csv"]

cur = mic[mic[CUR_MIC_VAL].notna() & mic["SMILES"].notna()]
ext = mic[mic[EXT_MIC_VAL].notna() & mic["SMILES"].notna()]

cur_by_sp = (
    cur.groupby(CUR_SPECIES)["SMILES"].nunique()
    .rename("n_mic_curated").reset_index()
    .rename(columns={CUR_SPECIES: "species"})
)
ext_by_sp = (
    ext.groupby(EXT_SPECIES)["SMILES"].nunique()
    .rename("n_mic_extracted").reset_index()
    .rename(columns={EXT_SPECIES: "species"})
)
ic50_valid = ic50_df[ic50_df[IC50_VAL].notna() & ic50_df["SMILES"].notna()]
ic50_by_sp = (
    ic50_valid.groupby(IC50_SPECIES)["SMILES"].nunique()
    .rename("n_ic50").reset_index()
    .rename(columns={IC50_SPECIES: "species"})
)
acc_valid = acc_df[acc_df[ACC_VAL].notna() & acc_df["SMILES"].notna()]
acc_by_sp = (
    acc_valid.groupby(ACC_SPECIES)["SMILES"].nunique()
    .rename("n_accumulation").reset_index()
    .rename(columns={ACC_SPECIES: "species"})
)

cov = (
    cur_by_sp
    .merge(ext_by_sp,  on="species", how="outer")
    .merge(ic50_by_sp, on="species", how="outer")
    .merge(acc_by_sp,  on="species", how="outer")
    .fillna(0)
)
for col in ["n_mic_curated", "n_mic_extracted", "n_ic50", "n_accumulation"]:
    cov[col] = cov[col].astype(int)
cov["n_total"]     = cov[["n_mic_curated", "n_mic_extracted", "n_ic50", "n_accumulation"]].sum(axis=1)
cov["patho_code"]  = cov["species"].map(organism_to_patho).fillna("")
cov["in_pipeline"] = cov["patho_code"].isin(pipeline_pathos)
cov = cov.sort_values("n_mic_curated", ascending=False).reset_index(drop=True)

print(f"\n  Total species across all assays : {len(cov):,}")
print(f"  Species with curated MIC data   : {(cov['n_mic_curated'] > 0).sum():,}")
print(f"  Species with extracted MIC data : {(cov['n_mic_extracted'] > 0).sum():,}")
print(f"  Species with curated IC50 data  : {(cov['n_ic50'] > 0).sum():,}")
print(f"  Species with accumulation data  : {(cov['n_accumulation'] > 0).sum():,}")
print(f"  Pipeline species present        : {cov['in_pipeline'].sum():,}")

print(f"\n  Top 30 species by curated MIC SMILES (* = pipeline):")
hdr = f"  {'Species':<47} {'MIC cur':>8} {'MIC ext':>8} {'IC50':>8} {'Accum':>7}"
print(hdr)
print("  " + "-" * (len(hdr) - 2))
for _, r in cov.head(30).iterrows():
    flag = "*" if r["in_pipeline"] else " "
    print(f"  {flag} {r['species']:<45} {r['n_mic_curated']:>8,} {r['n_mic_extracted']:>8,} "
          f"{r['n_ic50']:>8,} {r['n_accumulation']:>7,}")

cov.to_csv(os.path.join(outdir, "organism_coverage.csv"), index=False)
print(f"\nSaved organism_coverage.csv")


# ─── section 4: SPARK vs CO-ADD overlap ──────────────────────────────────────

print("\n" + "=" * 65)
print("SECTION 4: InChIKey overlap — SPARK curated MIC vs CO-ADD processed")
print("=" * 65)

smiles_info_path = os.path.join(processed_dir, "00_smiles_info.csv")
coadd_info = pd.read_csv(smiles_info_path)
coadd_keys = set(coadd_info["inchikey"].dropna().unique())
print(f"\n  CO-ADD processed compounds: {len(coadd_keys):,} InChIKeys")

print("  Computing InChIKeys for ALL SPARK curated MIC SMILES...")
all_cur_smiles = set(cur["SMILES"].dropna().unique())
print(f"  Unique SMILES to convert: {len(all_cur_smiles):,}")
all_spark_keys = smiles_set_to_inchikeys(all_cur_smiles)
print(f"  Valid SPARK InChIKeys   : {len(all_spark_keys):,}")

shared_all = all_spark_keys & coadd_keys
coadd_only = coadd_keys - all_spark_keys
spark_only = all_spark_keys - coadd_keys

print(f"\n  Overall overlap (SPARK curated MIC vs CO-ADD processed):")
print(f"    SPARK curated MIC total : {len(all_spark_keys):,}")
print(f"    CO-ADD processed total  : {len(coadd_keys):,}")
print(f"    Shared                  : {len(shared_all):,}  ({100*len(shared_all)/len(all_spark_keys):.1f}% of SPARK)")
print(f"    SPARK-only              : {len(spark_only):,}  ({100*len(spark_only)/len(all_spark_keys):.1f}% of SPARK)")
print(f"    CO-ADD-only             : {len(coadd_only):,}  ({100*len(coadd_only)/len(coadd_keys):.1f}% of CO-ADD)")

print(f"\n  Per-species breakdown (* = pipeline):")
hdr2 = f"  {'Species':<47} {'SPARK':>8} {'Shared':>8} {'SPARK-only':>11} {'% shared':>9}"
print(hdr2)
print("  " + "-" * (len(hdr2) - 2))

overlap_rows = []
for species in cur[CUR_SPECIES].dropna().unique():
    df_sp = cur[cur[CUR_SPECIES] == species]
    spark_keys_sp = smiles_set_to_inchikeys(df_sp["SMILES"].dropna().unique())
    if not spark_keys_sp:
        continue
    shared = spark_keys_sp & coadd_keys
    pct    = 100 * len(shared) / len(spark_keys_sp)
    overlap_rows.append({
        "species":     species,
        "patho_code":  organism_to_patho.get(species, ""),
        "in_pipeline": organism_to_patho.get(species, "") in pipeline_pathos,
        "n_spark":     len(spark_keys_sp),
        "n_shared":    len(shared),
        "n_spark_only": len(spark_keys_sp - coadd_keys),
        "pct_shared":  round(pct, 1),
    })

spark_vs_coadd = (
    pd.DataFrame(overlap_rows)
    .sort_values("n_spark", ascending=False)
    .reset_index(drop=True)
)

for _, r in spark_vs_coadd.iterrows():
    flag = "* " if r["in_pipeline"] else "  "
    print(f"  {flag}{r['species']:<45} {int(r['n_spark']):>8,} {int(r['n_shared']):>8,} "
          f"{int(r['n_spark_only']):>11,} {r['pct_shared']:>8.1f}%")

spark_vs_coadd.to_csv(os.path.join(outdir, "spark_vs_coadd.csv"), index=False)
print(f"\nSaved spark_vs_coadd.csv")


# ─── section 5: plots ─────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("SECTION 5: Plots")
print("=" * 65)

import stylia

# Format: print | Style: ersilia
stylia.set_format("print")
stylia.set_style("ersilia")

pal = stylia.CategoricalPalette("ersilia")
col_curated, col_extracted = pal.get(2)

nc = stylia.NamedColors()


# ── Plot A: top-25 species by curated MIC ────────────────────────────────────

def plot_mic_coverage(ax, df):
    ys = np.arange(len(df))
    ax.barh(ys, df["n_mic_curated"].values,   color=col_curated,   label="curated MIC")
    ax.barh(ys, df["n_mic_extracted"].values, color=col_extracted, label="extracted MIC",
            left=df["n_mic_curated"].values)
    labels = [
        (f"{r['species'][:42]} ◀" if r["in_pipeline"] else r["species"][:45])
        for _, r in df.iterrows()
    ]
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.legend()
    stylia.label(ax, xlabel="Unique SMILES", ylabel="",
                 title="SPARK MIC coverage — all species (top 25, ◀ = pipeline)")


plot_df_a = (
    cov[cov["n_mic_curated"] > 0]
    .sort_values("n_mic_curated", ascending=True)
    .tail(25)
)
fig, axs = stylia.create_figure(1, 1)
plot_mic_coverage(axs.next(), plot_df_a)
stylia.save_figure(os.path.join(outdir, "spark_mic_coverage_all.png"))
print("  Saved spark_mic_coverage_all.png")


# ── Plot B: SPARK vs CO-ADD overlap (top-20 by n_spark) ──────────────────────

def plot_overlap(ax, df):
    ys = np.arange(len(df))
    ax.barh(ys, df["n_shared"].values,     color=nc.plum, label="shared with CO-ADD")
    ax.barh(ys, df["n_spark_only"].values, color=nc.mint, label="SPARK-only",
            left=df["n_shared"].values)
    labels = [
        (f"{r['species'][:42]} ◀" if r["in_pipeline"] else r["species"][:45])
        for _, r in df.iterrows()
    ]
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.legend()
    stylia.label(ax, xlabel="Unique InChIKeys", ylabel="",
                 title="SPARK vs CO-ADD overlap — top 20 species (◀ = pipeline)")


plot_df_b = (
    spark_vs_coadd[spark_vs_coadd["n_spark"] > 0]
    .sort_values("n_spark", ascending=True)
    .tail(20)
)
fig, axs = stylia.create_figure(1, 1)
plot_overlap(axs.next(), plot_df_b)
stylia.save_figure(os.path.join(outdir, "spark_vs_coadd_overlap.png"))
print("  Saved spark_vs_coadd_overlap.png")

print("\nDone.")
