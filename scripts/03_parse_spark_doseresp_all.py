import os
import pandas as pd
import sys

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from utils import binarizer, merge_replicas

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, '..', 'data'))
processed_dir = os.path.abspath(os.path.join(root, '..', 'processed'))
config_dir = os.path.abspath(os.path.join(root, '..', 'config'))

pathogen_code = "abaumannii"
dc = pd.read_csv(os.path.join(config_dir, 'pathogens.csv'))
pathogen_name = dc[dc['pathogen_code'] == pathogen_code]['search_text'].values[0]

CUTOFFS_CONC = [1,10,40]

print("Starting with primary screening data")

df = pd.read_csv(os.path.join(data_dir, "spark", "preprocessed", "spark_mic.csv"), low_memory=False)

# Separate Curated & Transformed from Extracted and Uploaded and filter by pathogen
print(df.shape)
df1 = df[["smiles"]+[col for col in df.columns if "Curated & Transformed" in col]]
print(df1.shape)
df1 = df1.dropna(how="all", subset=[col for col in df.columns if "Curated & Transformed" in col])
print(df1.shape)

df2 = df[["smiles"]+[col for col in df.columns if "Extracted & Uploaded" in col]]
print(df2.shape)
df2 = df2.dropna(how="all", subset=[col for col in df.columns if "Extracted & Uploaded" in col])
print(df2.shape)

# Select pathogen
df1 = df1[df1["Curated & Transformed MIC Data: Species"]==pathogen_name]
df2 = df2[df2["Extracted & Uploaded MIC Data: Species"]==pathogen_name]
print("-----------------------------------")
print("Working on pathogen: ", pathogen_name, pathogen_code)
print("Curated: ", len(set(df1["smiles"])), "Extracted: ", len(set(df2["smiles"])))
print("SMILES in Curated not in Extracted: ", len(set(df1["smiles"])-set(df2["smiles"])))
print("SMILES in Extracted not in Curated: ", len(set(df2["smiles"])-set(df1["smiles"])))

#Extra SMILES in MIC file:
merged = pd.read_csv(os.path.join(processed_dir, "spark_individual", f"{pathogen_code}_spark_doseresponse_merged.csv"))
print("-----------------------------------")
print("SMILES in files not in MIC Curated: ", len(set(merged["smiles"])-set(df1["smiles"])))
print("SMILES in MIC Curated not in files: ", len(set(df1["smiles"])-set(merged["smiles"])))
print("SMILES in files not in MIC Extracted: ", len(set(merged["smiles"])-set(df2["smiles"])))
print("SMILES in MIC Extracted not in files: ", len(set(df2["smiles"])-set(merged["smiles"])))

# Separate Curated by Source and compare with original data
sources = df1["Curated & Transformed MIC Data: Source"].unique()
abau_strains = {"novartis": "NB48015", "co-add": "ATCC 19606", "merck": "IID876", "achaogen": "AABA1060"}
subset_sources = ["novartis", "co-add", "merck", "achaogen"]
for subset in subset_sources:
    subset_source = [s for s in sources if f"{subset}" in s.lower()]
    subset_mic = df1[df1["Curated & Transformed MIC Data: Source"].isin(subset_source)]
    print("-----------------------------------")
    print(f"{subset} comparison")
    values_col = 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
    strain_col = 'Curated & Transformed MIC Data: Strain'
    strain = abau_strains[subset]

    subset_mic = subset_mic[subset_mic[strain_col]==strain]

    subset_mic= subset_mic[["smiles", values_col]].copy()
    cutoff_cols = []
    for c in CUTOFFS_CONC:
            subset_mic[f"mic_{c}um"] = subset_mic[values_col].apply(lambda x: binarizer(x, c))
            cutoff_cols += [f"mic_{c}um"]

    rows_with_none = subset_mic[subset_mic.isnull().any(axis=1)]
    print("SMILES that cannot be binarised: ", len(rows_with_none))
    df_dicts = []
    for c in cutoff_cols:
        result = merge_replicas(subset_mic, c)
        df_dicts += [result]
    smiles_lists = [d["smiles"] for d in df_dicts]
    if not all(s == smiles_lists[0] for s in smiles_lists):
        raise ValueError("The SMILES differ between cut-offs, please check")
    subset_mic = pd.DataFrame({"smiles": smiles_lists[0]})
    for d in df_dicts:
        for key, values in d.items():
            if key != "smiles":
                subset_mic[key] = values
    print("Total datapoints: ", len(subset_mic))

    subset_orig = pd.read_csv(os.path.join(processed_dir, "spark_individual", f"{pathogen_code}_spark_doseresponse_{subset}.csv"))
    diff = pd.concat([subset_mic, subset_orig]).drop_duplicates(keep=False)
    if len(diff)==0:
        print(f"No differences between {subset} and MIC Curated Data")
    else:
        print(diff)

# If we remove the individual files data, do we have duplicates from outher sources with the individual files?
sources_names=[s for s in sources if any(f"{subset}" in s.lower() for subset in subset_sources)]
df1_no_ind = df1[~df1["Curated & Transformed MIC Data: Source"].isin(sources_names)]
print("-----------------------------------")
print("Remaining total datapoints without individual file sources: ", len(df1_no_ind))
print("Remaining unique datapoints without individual file sources: ", len(set(df1_no_ind["smiles"])))
print("Duplicates with the individual file sources: ", len(set(df1_no_ind["smiles"]).intersection(set(merged["smiles"]))))

# To keep data clean, we will not take into account the duplicates between the individual files and the other sources, and simply append the other sources

# Append other sources and merge without considering strain
df1_no_ind = df1[~df1["smiles"].isin(merged["smiles"])]
values_col = 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
strain_col = 'Curated & Transformed MIC Data: Strain'

df1_rest= df1_no_ind[["smiles", values_col]].copy()
cutoff_cols = []
for c in CUTOFFS_CONC:
        df1_rest[f"mic_{c}um"] = df1_rest[values_col].apply(lambda x: binarizer(x, c))
        cutoff_cols += [f"mic_{c}um"]
rows_with_none = df1_rest[df1_rest.isnull().any(axis=1)]
print("-----------------------------------")
print("Processing remaining CURATED data")
print("SMILES that cannot be binarised: ", len(rows_with_none))
df_dicts = []
for c in cutoff_cols:
    result = merge_replicas(df1_rest, c)
    df_dicts += [result]
smiles_lists = [d["smiles"] for d in df_dicts]
if not all(s == smiles_lists[0] for s in smiles_lists):
    raise ValueError("The SMILES differ between cut-offs, please check")
df1_rest = pd.DataFrame({"smiles": smiles_lists[0]})
for d in df_dicts:
    for key, values in d.items():
        if key != "smiles":
            df1_rest[key] = values
print("Total datapoints: ", len(df1_rest))
df1_rest.to_csv(os.path.join(processed_dir, "spark_individual", f"{pathogen_code}_spark_doseresponse_curated_rest.csv"), index=False)

# EXTRACT&UPLOAD
print("-----------------------------------")
print("Processing EXTRACTED data")

# Keep only smiles not already in individual files
df2_no_ind = df2[~df2["smiles"].isin(merged["smiles"])]
print("Datapoints in Extracted not present in individual files: ", df2_no_ind.shape)
print("Unique SMILES present in Extracted not in individual files: ",len(set(df2_no_ind["smiles"])))

# Compare results with CURATED file
df2_curated = df2_no_ind[df2_no_ind["smiles"].isin(df1_rest["smiles"])]
print("Datapoints in Extracted present in Curated Rest: ",df2_curated.shape)
print("Unique SMILES present in Extracted and in Curated Rest: ",len(set(df2_curated["smiles"])))

values_col = 'Extracted & Uploaded MIC Data: MIC (in µM) (µM)'
strain_col = 'Extracted & Uploaded MIC Data: Strain'

df2_rest= df2_curated[["smiles", values_col]].copy()
cutoff_cols = []
for c in CUTOFFS_CONC:
        df2_rest[f"mic_{c}um"] = df2_rest[values_col].apply(lambda x: binarizer(x, c))
        cutoff_cols += [f"mic_{c}um"]
rows_with_none = df2_rest[df2_rest.isnull().any(axis=1)]
print("-----------------------------------")
df_dicts = []
for c in cutoff_cols:
    result = merge_replicas(df2_rest, c)
    df_dicts += [result]
smiles_lists = [d["smiles"] for d in df_dicts]
if not all(s == smiles_lists[0] for s in smiles_lists):
    raise ValueError("The SMILES differ between cut-offs, please check")
df2_rest = pd.DataFrame({"smiles": smiles_lists[0]})
for d in df_dicts:
    for key, values in d.items():
        if key != "smiles":
            df2_rest[key] = values

df1_rest_in_extr = df1_rest[df1_rest["smiles"].isin(df2_rest["smiles"])]
diff = pd.concat([df1_rest_in_extr, df2_rest]).drop_duplicates(keep=False)
if len(diff)==0:
    print(f"No differences between Curated and Extracted MIC Data")
else:
    print("SMILES present in Curated and Extracted that differ: ", len(set(diff["smiles"])))

# Keep only SMILES not present in CURATED as the datasets are very similar
df2_rest = df2_no_ind[~df2_no_ind["smiles"].isin(df1_rest["smiles"])]
values_col = 'Extracted & Uploaded MIC Data: MIC (in µM) (µM)'
strain_col = 'Extracted & Uploaded MIC Data: Strain'
df2_rest= df2_rest[["smiles", values_col]].copy()
cutoff_cols = []
for c in CUTOFFS_CONC:
        df2_rest[f"mic_{c}um"] = df2_rest[values_col].apply(lambda x: binarizer(x, c))
        cutoff_cols += [f"mic_{c}um"]
rows_with_none = df2_rest[df2_rest.isnull().any(axis=1)]
print("-----------------------------------")
df_dicts = []
for c in cutoff_cols:
    result = merge_replicas(df2_rest, c)
    df_dicts += [result]
smiles_lists = [d["smiles"] for d in df_dicts]
if not all(s == smiles_lists[0] for s in smiles_lists):
    raise ValueError("The SMILES differ between cut-offs, please check")
df2_rest = pd.DataFrame({"smiles": smiles_lists[0]})
for d in df_dicts:
    for key, values in d.items():
        if key != "smiles":
            df2_rest[key] = values

print("Total datapoints added by Extracted: ", len(df2_rest))
df2_rest.to_csv(os.path.join(processed_dir, "spark_individual", f"{pathogen_code}_spark_doseresponse_extracted_rest.csv"), index=False)

# Join All Files
print("-----------------------------------")
print("Joining all files")
df_all = pd.concat([merged, df1_rest, df2_rest])
print("Total datapoints: ", len(df_all))
print("Total unique datapoints: ", len(set(df_all["smiles"])))

for c in CUTOFFS_CONC:
    df_ = df_all[["smiles", f"mic_{c}um"]]
    df_= df_[df_[f'mic_{c}um'].notna()]
    print(f"Final datapoints at cutoff {c} uM: ", len(df_))
    print("Positives: ", len(df_[df_[f'mic_{c}um']==1]), "Negatives: ", len(df_[df_[f'mic_{c}um']==0]))
    df_.to_csv(os.path.join(processed_dir, f"{pathogen_code}_spark_doseresp_{c}um.csv"), index=False)
