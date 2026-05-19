"""
CoADD analysis: cytotoxicity datasets, activity overlap plots, and dataset selection.

Part 1 — Cytotoxicity datasets
  Extract CC50 (HEK293 cells) and HC10 (Red blood cells) from the dose-response file.
  Binarise at [10, 25, 50] µM with direction = -1 (lower = more cytotoxic = active).
  Outputs:
    data/processed/coadd/06_citotox/cc50.csv
    data/processed/coadd/06_citotox/hc10.csv
  Columns: std_smiles, inchikey, mw, value, std, operator, replicas, cc50_10/25/50 (or hc10_*)

Part 2 — Bioactivity: inhibition % vs MIC per pathogen  ({patho}_bioactivity.png)
  For each pathogen that has both binarised inhibition and MIC files, find strains present
  in both. One panel per strain: inhibition value (x) vs MIC value (y), coloured by
  4-class activity (mid-cutoffs: inhib_75, mic_25). Spearman r annotation.

Part 3 — Cytotoxicity: MIC vs CC50 and inhibition % vs CC50 per pathogen
  ({patho}_cytotoxicity.png)
  Row 0: MIC (µM) vs CC50 (µM) per strain, coloured by activity (mic_25, cc50_25)
  Row 1: Inhibition % vs CC50 (µM) per strain, coloured by activity (inhib_75, cc50_25)
  Only rows with valid continuous values are plotted; Spearman r annotation.

Part 4 — Dataset selection
  Loads summary CSVs from scripts 03 and 05, flags datasets that should be excluded:
    small           : n_compounds < MIN_COMPOUNDS
    redundant_merged: merged where dominant strain covers >= COVERAGE_THRESHOLD of set
  Cutoff columns renamed to generic low/mid/high tiers from coadd_cutoffs.csv.
  Outputs:
    output/06_coadd_analysis/selection_table.csv
    output/06_coadd_analysis/{inhib,mic}_dataset_sizes.png
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(root, ".."))

from src.default import COADD_MIC_PATHO, COADD_INH_PATHO
from src.utils import extract_operator

# ---------------------------------------------------------------------------
# Parameters (dataset selection — Part 4)
# ---------------------------------------------------------------------------
MIN_COMPOUNDS      = 100   # datasets smaller than this are flagged
COVERAGE_THRESHOLD = 0.85  # merged flagged if dominant strain covers >= this fraction

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
config_dir  = os.path.join(root, "..", "config", "manual")
data_dir    = os.path.join(root, "..", "data", "raw", "coadd")
citotox_dir = os.path.join(root, "..", "data", "processed", "coadd", "06_citotox")
inh_bin_dir = os.path.join(root, "..", "data", "processed", "coadd", "03_binarised_inhibition")
mic_bin_dir = os.path.join(root, "..", "data", "processed", "coadd", "05_binarised_mic")
output_dir  = os.path.join(root, "..", "output", "06_coadd_analysis")

os.makedirs(citotox_dir, exist_ok=True)
os.makedirs(output_dir,  exist_ok=True)

# SMILES lookup
smiles_info = pd.read_csv(
    os.path.join(root, "..", "data", "processed", "coadd", "00_smiles_info.csv")
)
smiles_lookup = smiles_info.set_index("smiles")[["std_smiles", "inchikey", "mw"]].to_dict("index")

# ---------------------------------------------------------------------------
# Part 1 — Cytotoxicity datasets
# ---------------------------------------------------------------------------
CITOTOX_CUTOFFS  = [10, 25, 50]
CC50_CUTOFF_COLS = [f"cc50_{c}" for c in CITOTOX_CUTOFFS]
HC10_CUTOFF_COLS = [f"hc10_{c}" for c in CITOTOX_CUTOFFS]


def binarize_citotox(operator, value, cutoff):
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


def binarise_citotox(df, cutoffs, cutoff_cols):
    if df.empty:
        return pd.DataFrame(columns=["std_smiles", "inchikey", "mw",
                                     "value", "std", "operator", "replicas"] + cutoff_cols)
    rows = []
    for std_smiles, grp in df.groupby("std_smiles"):
        first   = grp.iloc[0]
        vals    = grp["value"].dropna().tolist()
        n       = len(vals)
        avg_val = round(float(np.mean(vals)), 3) if vals else np.nan
        std_val = round(float(np.std(vals, ddof=1)), 3) if n > 1 else 0.0
        ops     = grp["operator"].dropna().unique().tolist()
        agg_op  = ops[0] if len(ops) == 1 else "mixed"
        row = {
            "std_smiles": std_smiles,
            "inchikey":   first["inchikey"],
            "mw":         round(float(first["mw"]), 3) if not pd.isna(first["mw"]) else np.nan,
            "value":      avg_val,
            "std":        std_val,
            "operator":   agg_op,
            "replicas":   n,
        }
        for cutoff, col in zip(cutoffs, cutoff_cols):
            labels = [
                binarize_citotox(r["operator"], r["value"], cutoff)
                for _, r in grp.iterrows()
                if not pd.isna(r["value"])
            ]
            row[col] = aggregate_labels(labels) if labels else -1
        rows.append(row)
    return pd.DataFrame(rows)[
        ["std_smiles", "inchikey", "mw", "value", "std", "operator", "replicas"] + cutoff_cols
    ].reset_index(drop=True)


def mic_to_um(numeric, unit, mw):
    if pd.isna(numeric) or pd.isna(mw) or mw <= 0:
        return None
    if unit == "uM":
        return float(numeric)
    if unit == "ug/mL":
        return (float(numeric) / mw) * 1000
    return None


print("=" * 60)
print("PART 1 — Cytotoxicity datasets")
print("=" * 60)

dr = pd.read_csv(
    os.path.join(data_dir, "CO-ADD_DoseResponseData_r03_01-02-2020_CSV.csv"),
    low_memory=False,
)
hs = dr[
    (dr["ORGANISM"] == "Homo sapiens") &
    dr["SMILES"].notnull() &
    dr["DRVAL_MEDIAN"].notnull()
].copy()

for drval_type, cutoff_cols in [("CC50", CC50_CUTOFF_COLS), ("HC10", HC10_CUTOFF_COLS)]:
    df_t = hs[hs["DRVAL_TYPE"] == drval_type].copy()
    print(f"\n[{drval_type}] {len(df_t):,} rows | units: {df_t['DRVAL_UNIT'].value_counts().to_dict()}")

    df_t["std_smiles"] = df_t["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("std_smiles"))
    df_t["inchikey"]   = df_t["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("inchikey"))
    df_t["mw"]         = df_t["SMILES"].map(lambda s: (smiles_lookup.get(s) or {}).get("mw"))

    parsed           = df_t["DRVAL_MEDIAN"].apply(lambda v: extract_operator(str(v)))
    df_t["operator"] = parsed.apply(lambda t: t[0])
    df_t["numeric"]  = parsed.apply(lambda t: float(t[1]) if t[1] is not None else None)
    df_t["value"]    = df_t.apply(
        lambda r: mic_to_um(r["numeric"], r["DRVAL_UNIT"], r["mw"]), axis=1
    )
    df_t = df_t.dropna(subset=["std_smiles", "value"]).reset_index(drop=True)

    df_bin = binarise_citotox(df_t, CITOTOX_CUTOFFS, cutoff_cols)
    fname  = f"{drval_type.lower()}.csv"
    df_bin.to_csv(os.path.join(citotox_dir, fname), index=False)

    n = len(df_bin)
    for col in cutoff_cols:
        definitive = df_bin[df_bin[col] != -1]
        n_act = int((definitive[col] == 1).sum())
        n_def = len(definitive)
        rate  = n_act / n_def if n_def else float("nan")
        n_inc = int((df_bin[col] == -1).sum())
        print(f"  {col}: n={n}, active={n_act}, inconclusive={n_inc}, rate={rate:.1%}")
    print(f"  -> saved {fname}")

cc50 = pd.read_csv(os.path.join(citotox_dir, "cc50.csv"))
hc10 = pd.read_csv(os.path.join(citotox_dir, "hc10.csv"))

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
import stylia

# Format: print | Style: ersilia
stylia.set_format("print")
stylia.set_style("ersilia")

nc = stylia.NamedColors()

# Category colors: (both, x_only, y_only, neither)
CAT_COLORS = [nc.plum, nc.orange, nc.mint, nc.gray]


def figure_width(n_cols):
    return 0.5 if n_cols == 1 else 1.0


def assign_colors(df, x_act_col, y_act_col):
    colors = []
    for _, row in df.iterrows():
        xa = int(row[x_act_col]) == 1
        ya = int(row[y_act_col]) == 1
        if xa and ya:
            colors.append(CAT_COLORS[0])
        elif xa:
            colors.append(CAT_COLORS[1])
        elif ya:
            colors.append(CAT_COLORS[2])
        else:
            colors.append(CAT_COLORS[3])
    return colors


def draw_scatter(ax, x_vals, y_vals, colors, cat_labels, xlabel, ylabel, title):
    """Scatter with Spearman r annotation and 4-class legend."""
    valid = ~(np.isnan(x_vals) | np.isnan(y_vals))
    if valid.sum() >= 3:
        x_v = x_vals[valid]
        y_v = y_vals[valid]
        c_v = [colors[i] for i, v in enumerate(valid) if v]
        ax.scatter(x_v, y_v, color=c_v)
        r, _ = spearmanr(x_v, y_v)
        ax.annotate(f"r = {r:.2f}", xy=(0.05, 0.92), xycoords="axes fraction")
    # Legend
    for label, color in cat_labels.items():
        ax.scatter([], [], color=color, label=label)
    ax.legend()
    stylia.label(ax, xlabel=xlabel, ylabel=ylabel, title=title)


# ---------------------------------------------------------------------------
# Part 2 — Bioactivity: inhibition % vs MIC per pathogen
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("PART 2 — Bioactivity (Inh % vs MIC)")
print("=" * 60)

INH_MID_COL = "inhib_75"
MIC_MID_COL = "mic_25"

INH_MIC_CATS = {
    "Both active":  CAT_COLORS[0],
    "Inh active":   CAT_COLORS[1],
    "MIC active":   CAT_COLORS[2],
    "Neither":      CAT_COLORS[3],
}

for patho_code in sorted(set(COADD_INH_PATHO) & set(COADD_MIC_PATHO)):
    inh_files = {
        f.replace(f"{patho_code}_", "").replace(".csv", ""): f
        for f in os.listdir(inh_bin_dir)
        if f.startswith(f"{patho_code}_") and not f.endswith("merged.csv")
    }
    mic_files = {
        f.replace(f"{patho_code}_", "").replace(".csv", ""): f
        for f in os.listdir(mic_bin_dir)
        if f.startswith(f"{patho_code}_") and not f.endswith("merged.csv")
    }
    common_strains = sorted(set(inh_files) & set(mic_files))
    if not common_strains:
        print(f"[{patho_code}] No common strains, skipping.")
        continue

    print(f"\n[{patho_code}] {len(common_strains)} strain(s): {common_strains}")

    n_cols = len(common_strains)
    fig, axs = stylia.create_figure(1, n_cols, width=figure_width(n_cols))

    for strain in common_strains:
        df_inh = pd.read_csv(os.path.join(inh_bin_dir, inh_files[strain]))
        df_mic = pd.read_csv(os.path.join(mic_bin_dir, mic_files[strain]))

        merged = (
            df_inh[["std_smiles", "value", INH_MID_COL]]
            .rename(columns={"value": "inh_value"})
            .merge(
                df_mic[["std_smiles", "value", MIC_MID_COL]]
                .rename(columns={"value": "mic_value"}),
                on="std_smiles",
            )
            .dropna(subset=["inh_value", "mic_value"])
        )
        print(f"  {strain}: {len(merged)} overlapping molecules")

        colors = assign_colors(merged, INH_MID_COL, MIC_MID_COL)
        ax = axs.next()
        draw_scatter(
            ax,
            merged["inh_value"].values.astype(float),
            merged["mic_value"].values.astype(float),
            colors, INH_MIC_CATS,
            xlabel="Inhibition (%)",
            ylabel="MIC (µM)",
            title=f"{patho_code} / {strain}",
        )

    stylia.save_figure(os.path.join(output_dir, f"{patho_code}_bioactivity.png"))
    plt.close("all")
    print(f"  -> {patho_code}_bioactivity.png")

# ---------------------------------------------------------------------------
# Part 3 — Cytotoxicity: MIC vs CC50 and Inh % vs CC50 per pathogen
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("PART 3 — Cytotoxicity (MIC vs CC50 | Inh % vs CC50)")
print("=" * 60)

CC50_MID_COL = "cc50_25"

MIC_CC50_CATS = {
    "Both active":   CAT_COLORS[0],
    "MIC active":    CAT_COLORS[1],
    "CC50 active":   CAT_COLORS[2],
    "Neither":       CAT_COLORS[3],
}
INH_CC50_CATS = {
    "Both active":   CAT_COLORS[0],
    "Inh active":    CAT_COLORS[1],
    "CC50 active":   CAT_COLORS[2],
    "Neither":       CAT_COLORS[3],
}

all_pathogens = sorted(set(COADD_MIC_PATHO + COADD_INH_PATHO))

for patho_code in all_pathogens:
    # Row 0: MIC vs CC50
    mic_cc50_panels = []
    for f in sorted(os.listdir(mic_bin_dir)):
        if not f.startswith(f"{patho_code}_") or f.endswith("merged.csv"):
            continue
        strain = f.replace(f"{patho_code}_", "").replace(".csv", "")
        df_mic = pd.read_csv(os.path.join(mic_bin_dir, f))
        merged = (
            df_mic[["std_smiles", "value", MIC_MID_COL]]
            .rename(columns={"value": "mic_value"})
            .merge(
                cc50[["std_smiles", "value", CC50_MID_COL]]
                .rename(columns={"value": "cc50_value"}),
                on="std_smiles",
            )
            .dropna(subset=["mic_value", "cc50_value"])
        )
        if len(merged) >= 3:
            mic_cc50_panels.append((strain, merged))

    # Row 1: Inh % vs CC50
    inh_cc50_panels = []
    for f in sorted(os.listdir(inh_bin_dir)):
        if not f.startswith(f"{patho_code}_") or f.endswith("merged.csv"):
            continue
        strain = f.replace(f"{patho_code}_", "").replace(".csv", "")
        df_inh = pd.read_csv(os.path.join(inh_bin_dir, f))
        merged = (
            df_inh[["std_smiles", "value", INH_MID_COL]]
            .rename(columns={"value": "inh_value"})
            .merge(
                cc50[["std_smiles", "value", CC50_MID_COL]]
                .rename(columns={"value": "cc50_value"}),
                on="std_smiles",
            )
            .dropna(subset=["inh_value", "cc50_value"])
        )
        if len(merged) >= 3:
            inh_cc50_panels.append((strain, merged))

    if not mic_cc50_panels and not inh_cc50_panels:
        print(f"[{patho_code}] No CC50 overlap, skipping.")
        continue

    n_rows = (1 if not mic_cc50_panels else 1) + (1 if inh_cc50_panels else 0)
    n_cols = max(len(mic_cc50_panels), len(inh_cc50_panels), 1)

    print(f"\n[{patho_code}] {len(mic_cc50_panels)} MIC-CC50 panel(s), {len(inh_cc50_panels)} Inh-CC50 panel(s)")

    fig, axs = stylia.create_figure(n_rows, n_cols, width=figure_width(n_cols))

    # Row 0 — MIC vs CC50
    for strain, merged in mic_cc50_panels:
        print(f"  MIC-CC50 {strain}: {len(merged)} molecules")
        colors = assign_colors(merged, MIC_MID_COL, CC50_MID_COL)
        ax = axs.next()
        draw_scatter(
            ax,
            merged["mic_value"].values.astype(float),
            merged["cc50_value"].values.astype(float),
            colors, MIC_CC50_CATS,
            xlabel="MIC (µM)",
            ylabel="CC50 (µM)",
            title=f"{patho_code} / {strain}",
        )
    # Blank remaining cells in row 0
    for _ in range(n_cols - len(mic_cc50_panels)):
        ax = axs.next()
        ax.axis("off")

    # Row 1 — Inh % vs CC50 (only if present)
    if inh_cc50_panels:
        for strain, merged in inh_cc50_panels:
            print(f"  Inh-CC50 {strain}: {len(merged)} molecules")
            colors = assign_colors(merged, INH_MID_COL, CC50_MID_COL)
            ax = axs.next()
            draw_scatter(
                ax,
                merged["inh_value"].values.astype(float),
                merged["cc50_value"].values.astype(float),
                colors, INH_CC50_CATS,
                xlabel="Inhibition (%)",
                ylabel="CC50 (µM)",
                title=f"{patho_code} / {strain}",
            )
        for _ in range(n_cols - len(inh_cc50_panels)):
            ax = axs.next()
            ax.axis("off")

    stylia.save_figure(os.path.join(output_dir, f"{patho_code}_cytotoxicity.png"))
    plt.close("all")
    print(f"  -> {patho_code}_cytotoxicity.png")

# ---------------------------------------------------------------------------
# Part 4 — Dataset selection
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("PART 4 — Dataset selection")
print("=" * 60)

cutoffs_cfg = pd.read_csv(os.path.join(config_dir, "coadd_cutoffs.csv"))


def build_rename_map(assay_type, prefix):
    """Map summary CSV column names → generic low/mid/high names."""
    row = cutoffs_cfg[cutoffs_cfg["assay_type"] == assay_type].iloc[0]
    tiers = {
        "low":  int(row["cutoff_low"]),
        "mid":  int(row["cutoff_mid"]),
        "high": int(row["cutoff_high"]),
    }
    return {
        f"n_active_{prefix}_{val}":    f"n_active_{tier}"
        for tier, val in tiers.items()
    } | {
        f"active_rate_{prefix}_{val}": f"active_rate_{tier}"
        for tier, val in tiers.items()
    }


INH_RENAME = build_rename_map("inhib", "inhib")
MIC_RENAME = build_rename_map("mic",   "mic")

inh_sum = pd.read_csv(
    os.path.join(root, "..", "output", "03_binarise_coadd_inhibition", "summary.csv")
).rename(columns=INH_RENAME)
mic_sum = pd.read_csv(
    os.path.join(root, "..", "output", "05_binarise_coadd_mic", "summary.csv")
).rename(columns=MIC_RENAME)

inh_sum["assay_type"] = "inhib"
mic_sum["assay_type"] = "mic"

CUTOFF_COLS = ["n_active_low", "n_active_mid", "n_active_high",
               "active_rate_low", "active_rate_mid", "active_rate_high"]
base_cols   = ["patho_code", "strain_code", "assay_type", "n_compounds"]

sel = pd.concat(
    [inh_sum[base_cols + CUTOFF_COLS], mic_sum[base_cols + CUTOFF_COLS]],
    ignore_index=True,
)

# Coverage: for merged rows, fraction of merged covered by the dominant strain
coverage_map = {}
for (patho, assay), grp in sel.groupby(["patho_code", "assay_type"]):
    merged_row = grp[grp["strain_code"] == "merged"]
    if merged_row.empty:
        continue
    n_merged   = int(merged_row["n_compounds"].values[0])
    n_dominant = int(grp[grp["strain_code"] != "merged"]["n_compounds"].max())
    if n_merged > 0:
        coverage_map[(patho, assay)] = round(n_dominant / n_merged, 4)

sel["_coverage"] = sel.apply(
    lambda r: coverage_map.get((r["patho_code"], r["assay_type"]), np.nan)
    if r["strain_code"] == "merged" else np.nan,
    axis=1,
)


def assign_flag(row):
    if row["n_compounds"] < MIN_COMPOUNDS:
        return "small"
    if row["strain_code"] == "merged" and not np.isnan(row["_coverage"]) \
            and row["_coverage"] >= COVERAGE_THRESHOLD:
        return "redundant_merged"
    return "ok"


sel["flag"] = sel.apply(assign_flag, axis=1)
sel["keep"] = sel["flag"] == "ok"

print(f"\nParameters: MIN_COMPOUNDS={MIN_COMPOUNDS}, COVERAGE_THRESHOLD={COVERAGE_THRESHOLD}")
print(f"Keep: {sel['keep'].sum()} / {len(sel)} datasets\n")
print("Excluded:")
print(sel[~sel["keep"]][["patho_code", "assay_type", "strain_code",
                          "n_compounds", "flag"]].to_string(index=False))

col_order = base_cols + CUTOFF_COLS + ["flag", "keep"]
sel[col_order].to_csv(os.path.join(output_dir, "selection_table.csv"), index=False)
print(f"\n-> selection_table.csv")

# Bar charts — one per assay type
FLAG_COLOR = {"ok": nc.blue, "small": nc.orange, "redundant_merged": nc.gray}
FLAG_LABEL = {
    "ok":               "Selected",
    "small":            f"Too small (< {MIN_COMPOUNDS})",
    "redundant_merged": f"Redundant merged (dominant >= {COVERAGE_THRESHOLD:.0%})",
}

for assay_type in ["inhib", "mic"]:
    sub = sel[sel["assay_type"] == assay_type].copy()
    sub["_is_merged"] = (sub["strain_code"] == "merged").astype(int)
    sub = sub.sort_values(["patho_code", "_is_merged", "strain_code"]).reset_index(drop=True)

    labels = (sub["patho_code"] + " / " + sub["strain_code"]).tolist()
    fig, axs = stylia.create_figure(1, 1)
    ax = axs.next()
    y_pos = np.arange(len(sub))
    ax.barh(y_pos, sub["n_compounds"].tolist(), color=[FLAG_COLOR[f] for f in sub["flag"]])
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.axvline(MIN_COMPOUNDS, linestyle="--", color=nc.gray)
    handles = [
        Patch(facecolor=FLAG_COLOR[flag], label=FLAG_LABEL[flag])
        for flag in ["ok", "small", "redundant_merged"]
        if flag in sub["flag"].values
    ]
    ax.legend(handles=handles)
    stylia.label(ax, xlabel="Compounds", ylabel="", title=f"{assay_type} datasets")
    stylia.save_figure(os.path.join(output_dir, f"{assay_type}_dataset_sizes.png"))
    plt.close("all")
    print(f"-> {assay_type}_dataset_sizes.png")

print("\nDone.")
