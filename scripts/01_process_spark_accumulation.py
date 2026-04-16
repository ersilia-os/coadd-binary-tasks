import os
import sys
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)

datapath = os.path.join(root,"..","data", "spark","preprocessed")
resultspath = os.path.join(root,"..","processed", "spark","accumulation")
if not os.path.exists(resultspath):
    os.makedirs(resultspath)

df = pd.read_csv(os.path.join(datapath, "spark_accumulation.csv"))

print(df.shape)

exp_cols = ['Curated & Transformed Accumulation Data: Accumulated Compound (nmol/1E12 CFU)', 'Curated & Transformed Accumulation Data: Accumulated Compound (pg/1E9 CFU)',
            'Curated & Transformed Accumulation Data: Accumulated Compound (ng/mL)','Curated & Transformed Accumulation Data: Accumulated Compound (μM)']

df = df.dropna(subset=exp_cols, how="all")
print(df.shape)
duplicates = df[df["Compound Name"].duplicated()]["Compound Name"].unique()
print("Keeping only rows with data, we have: ", len(duplicates), " duplicated smiles")

# Check publications
publ_col = "Curated & Transformed Accumulation Data: DOI"
publications = set(df[publ_col].tolist())
print(publications)
print("Data not associated to a publication:", len(df[df[publ_col].isna()]))



"""
pheno = list(set(df["Curated & Transformed Accumulation Data: Accumulation phenotype"].tolist()))
pheno_discard = ['Efflux overexpressor', 'Hyperpermeable', 'Efflux deficient', 'Efflux deficient; Hyperpermeable']

# From Richter et al, 2017, get the controls
poscontrol = ["tetracycline", "ciprofloxacin", "chloramphenicol"]
negcontrol = ["novobiocin", "erythromycin", "rifampicin", "vancomycin", "daptomycin", "clindamycin", "mupirocin", "fusidic acid", "ampicillin"]

cols_to_keep = {
    'smiles':'smiles',
    'Compound Name': "id", 
    'Synonyms': "synonym", 
    'Curated & Transformed Accumulation Data: Accumulated Compound (nmol/1E12 CFU)': 'acc_nmol/1E12CFU',
    'Curated & Transformed Accumulation Data: Accumulated Compound (μM)': 'acc_um',
    'Curated & Transformed Accumulation Data: Test article type': "type",
    'Curated & Transformed Accumulation Data: Compound incubation concentration (μM)': "inc_um",
    'Curated & Transformed Accumulation Data: Assay incubation time (min)': "inc_min", 
    'Curated & Transformed Accumulation Data: Species':"species",
    'Curated & Transformed Accumulation Data: Strain':"strain",
}

# Separate by pathogen to be able to merge duplicates if it is the case
unique_species = df['Curated & Transformed Accumulation Data: Species'].dropna().unique().tolist()
for s in unique_species:
    print("Species: ", s)
    df_ = df[df['Curated & Transformed Accumulation Data: Species']==s]
    duplicates = df_[df_["Compound Name"].duplicated()]["Compound Name"].unique()
    print("Total Compounds:", len(df_), "Duplicates:", len(duplicates))

    # we do not want experiments that have extra additives
    df_ = df_[~df_["Curated & Transformed Accumulation Data: Additives"].notna()]
    duplicates = df_[df_["Compound Name"].duplicated()]["Compound Name"].unique()
    print("Keeping only experiments without additives, we have: ", len(duplicates), " duplicated smiles")

    # We are only focused on WT strains
    df_ = df_[~df_["Curated & Transformed Accumulation Data: Accumulation phenotype"].isin(pheno_discard)]
    duplicates = df_[df_["Compound Name"].duplicated()]["Compound Name"].unique()
    print("Keeping only WT phenotypes, we have: ", len(duplicates), " duplicated smiles")

    test_col = "Curated & Transformed Accumulation Data: Test article type"
    control_mask = df_[test_col].isin(["Control / Reference compound"])
    dups = df_[control_mask]["Compound Name"].duplicated().unique()
    to_drop = df_[control_mask].loc[
        df_[control_mask]["Compound Name"].duplicated(keep="first")
    ].index
    df_ = df_.drop(index=to_drop)
    duplicates = df_[df_["Compound Name"].duplicated()]["Compound Name"].unique()
    print("Remaining duplicates after removing controls/reference:", len(duplicates))
    #manually clean up remaining dups:
    if s == "Pseudomonas aeruginosa":
        compounds_to_fix = ["SPK-0151060", "SPK-0151061", "SPK-0151062"]
        mask_comp = df_["Compound Name"].isin(compounds_to_fix)
        to_drop = df_[mask_comp & (df_["Curated & Transformed Accumulation Data: Assay incubation time (min)"] != 30)].index
        df_ = df_.drop(index=to_drop)
        duplicates = df_[df_["Compound Name"].duplicated()]["Compound Name"].unique()
        print("Remaining duplicates after final cleaning:", len(duplicates))
    if s == "Escherichia coli":
        keep_time_ten = ["SPK-0005037", "SPK-0151060", "SPK-0151061", "SPK-0151062"]
        mask_comp = df_["Compound Name"].isin(keep_time_ten)
        to_drop = df_[mask_comp & (df_["Curated & Transformed Accumulation Data: Assay incubation time (min)"] != 10)].index
        df_ = df_.drop(index=to_drop)
        
        keep_time_fifteen = ["SPK-0244706"]
        mask_comp = df_["Compound Name"].isin(keep_time_fifteen)
        to_drop = df_[mask_comp & (df_["Curated & Transformed Accumulation Data: Assay incubation time (min)"] != 15)].index
        df_ = df_.drop(index=to_drop)
        
        keep_concentration_four = ["SPK-0255033", "SPK-0255037"]
        mask_comp = df_["Compound Name"].isin(keep_concentration_four)
        to_drop = df_[mask_comp & (df_["Curated & Transformed Accumulation Data: Compound incubation concentration (μg/mL)"] != 4)].index
        df_ = df_.drop(index=to_drop)

        average_exp = ["SPK-0004789", "SPK-0004808", "SPK-0004809", "SPK-0112980", "SPK-0131058", "SPK-0131057", "SPK-0255051", "SPK-0255058", "SPK-0255060",
                       "SPK-0005037", "SPK-0151060", "SPK-0151061", "SPK-0151062", "SPK-0255037"] 

        mask_avg = df_["Compound Name"].isin(average_exp)
        df_avg = df_[mask_avg].sort_values("Compound Name").copy()

        rows = []
        for name, g in df_avg.groupby("Compound Name"):
            base = g.iloc[0].copy()
            for col in exp_cols:
                base[col] = pd.to_numeric(g[col], errors="coerce").mean()
            rows.append(base)

        df_avg_single = pd.DataFrame(rows)
        df_ = df_[~mask_avg].copy()
        df_ = pd.concat([df_, df_avg_single], ignore_index=True)
        duplicates = df_[df_["Compound Name"].duplicated()]["Compound Name"].unique()
        print("Remaining duplicates after final cleaning:", len(duplicates))
    print("Final number of compounds: ", len(df_))

    keep_cols = []
    for k,v in cols_to_keep.items():
        keep_cols += [k]
    df_ = df_[keep_cols]
    df_ = df_.rename(columns = cols_to_keep)

# Binarization will be done based on Richter et al, 2017 (Fig 1A); threshold 300 nmol/CFA
# If in known lists, classify as 0 or 1 without looking at experimental data

    df_.to_csv(os.path.join(resultspath, f"accumulation_{s}.csv"), index=False)
"""