import numpy as np
import collections
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
from standardiser import standardise
import pandas as pd
import os

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, '..', 'data'))
processed_dir = os.path.abspath(os.path.join(root, '..', 'processed'))
config_dir = os.path.abspath(os.path.join(root, '..', 'config'))

pathogen_code = "abaumannii"
dc = pd.read_csv(os.path.join(config_dir, 'pathogens.csv'))
pathogen_name = dc[dc['pathogen_code'] == pathogen_code]['search_text'].values[0]

def standardise_smiles(df):
    std_smiles = []
    for smiles in tqdm(list(df["smiles"])):
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = standardise.run(mol)
            std_smiles.append(Chem.MolToSmiles(mol))
        except:
            std_smiles.append(None)
    df["std_smiles"] = std_smiles
    df = df[df["std_smiles"].notnull()]
    return df

CUTOFFS_PERC = [50, 80, 90]
CUTOFFS_CONC = [1,10,40]

print("Working on pathogen: ", pathogen_name, pathogen_code)

print("Starting with primary screening data")

# Read in the spark data (MIC)


mic_files = {
    "mk": [pd.read_csv(os.path.join(data_dir, "spark","SPARK Data Merck & Kyorin Contribution.csv"), low_memory=False),
           'Curated & Transformed MIC Data: MIC (in µM) (µM)',
           'Curated & Transformed MIC Data: Species'
           ],
    "qe" :[pd.read_csv(os.path.join(data_dir, "spark","SPARK Data Quave Lab {Emory University} Publications.csv"), low_memory=False),
            'Porras 2020 Chem. Rev. XX:XXXX: MIC (in µM) (µM)',
            'Porras 2020 Chem. Rev. XX:XXXX: Species'
            ],
    "no": [pd.read_csv(os.path.join(data_dir, "spark","SPARK Data Novartis Contribution.csv"), low_memory=False),
           'Curated & Transformed MIC Data: MIC (in µM) (µM)',
           'Curated & Transformed MIC Data: Species'
           ],
    "ac": [pd.read_csv(os.path.join(data_dir, "spark","SPARK Data Achaogen Contribution.csv"), low_memory=False),
              'Curated & Transformed MIC Data: MIC (in µM) (µM)',
                'Curated & Transformed MIC Data: Species'
              ],
    "co": [pd.read_csv(os.path.join(data_dir, "spark","SPARK Data CO-ADD Contribution.csv"), low_memory=False),
              'Curated & Transformed MIC Data: MIC (in µM) (µM)',
              'Curated & Transformed MIC Data: Species'
              ]
           }

mic = pd.read_csv(os.path.join(data_dir, "spark","SPARK MIC Data.csv"), low_memory=False)

# STEP 1: Treat mic files
for k,v in mic_files.items():
    v[0] = v[0][v[0]["SMILES"].notnull()]
    v[0] = v[0][v[0][v[1]].notnull()]
    v[0] = v[0][v[0][v[2]]==pathogen_name]
    v[0] = v[0][["SMILES", v[1]]]
    v[0] = v[0].rename(columns={"SMILES":"smiles", v[1]:"mic_um"})
    mic_files[k][0] = v[0]
    print(k, v[0].shape)

#STEP 2: Treat the main MIC file
mic = mic[mic["SMILES"].notnull()]
mic = mic[["SMILES", "Curated & Transformed MIC Data: Species", "Curated & Transformed MIC Data: MIC (in µM) (µM)", "Extracted & Uploaded MIC Data: Species", "Extracted & Uploaded MIC Data: MIC (in µM) (µM)"]]
mic = mic.rename(columns={"SMILES": "smiles", 
                          "Curated & Transformed MIC Data: Species": "curated_species",
                          "Curated & Transformed MIC Data: MIC (in µM) (µM)": "curated_mic_um",
                          "Extracted & Uploaded MIC Data: Species": "extracted_species",
                          "Extracted & Uploaded MIC Data: MIC (in µM) (µM)": "extracted_mic_um"})

print(mic.columns)
mic =mic[
    (mic["curated_species"] == pathogen_name) |
    (mic["extracted_species"] == pathogen_name)
]
mic = mic.drop(columns=[
    "curated_species",
    "extracted_species"
])
mic = mic[
    mic["curated_mic_um"].notna() |
    mic["extracted_mic_um"].notna()
]
mic["mic_match"] = (
    mic["curated_mic_um"].isna() |
    mic["extracted_mic_um"].isna() |
    (mic["curated_mic_um"] == mic["extracted_mic_um"]))

print(mic["mic_match"].value_counts())
print(" SMILES not Matching between CURATED and EXTRACTED that will be eliminated: ",len(mic[mic["mic_match"] == False]))
false_match_rows = mic[mic["mic_match"] == False]
print(false_match_rows.head(10))
mic = mic.assign(mic_um=None) 
mic["mic_um"] = mic.apply(
    lambda row: row["curated_mic_um"] if row["mic_match"] else None, axis=1
)
mic = mic[mic["mic_um"].notna()]
mic = mic.reset_index(drop=True)
mic = mic.drop(columns=["curated_mic_um", "extracted_mic_um", "mic_match"])
mic["mic_um"]=mic["mic_um"].astype(str)


# STEP 3: Compare the mic files with the main mic file
for k, v in mic_files.items():
    mic_file_df = v[0]
    v[0]["mic_um"] =  v[0]["mic_um"].astype(str)
    matched_rows = mic_file_df.merge(mic, on=["smiles", "mic_um"], how="inner")
    num_matches = len(matched_rows)
    num_mismatches = len(mic_file_df) - num_matches
    print(f"File {k}:")
    print(f"  Total rows: {len(mic_file_df)}")
    print(f"  Matches in mic: {num_matches}")
    print(f"  Not in mic: {num_mismatches}\n")

combined_mic_files = pd.concat([v[0] for v in mic_files.values()], ignore_index=True)
combined_mic_files["smiles"] = combined_mic_files["smiles"].astype(str)
combined_mic_files["mic_um"] = combined_mic_files["mic_um"].astype(str)
mic["smiles"] = mic["smiles"].astype(str)
mic["mic_um"] = mic["mic_um"].astype(str)
matching_rows = mic.merge(combined_mic_files, on=["smiles", "mic_um"], how="inner")
different_rows = mic.merge(combined_mic_files, on=["smiles", "mic_um"], how="outer", indicator=True)
different_rows = different_rows[different_rows["_merge"] != "both"]
print(f"Total rows in mic: {len(mic)}")
print(f"Total rows in combined mic_files: {len(combined_mic_files)}")
print(f"Matching rows: {len(matching_rows)}")
print(f"Different rows: {len(different_rows)}")

# STEP 4: Concatenate dataframes and drop exact duplicates
all_mic = pd.concat([mic, combined_mic_files], ignore_index=True)
all_mic = all_mic.drop_duplicates()

# STEP 5: Standardise SMILES

all_mic = standardise_smiles(all_mic)
all_mic = all_mic[["std_smiles", "mic_um"]]
all_mic.rename(columns={"std_smiles": "smiles"}, inplace=True)

# STEP 6: Convert units to binary
def extract_operator(v):
    v = str(v)
    try:
        float(v)
        return "="
    except:
        if v.startswith(">"):
            return ">"
        elif v.startswith("<"):
            return "<"
        return None

def binarizer(v, cutoff):
    operator = extract_operator(v)
    if operator is None:
        return None
    if operator == "=":
        v = float(v)
        if v <= cutoff:
            return 1
        else:
            return 0
    if operator == ">":
        v = "".join([x for x in str(v) if x.isdigit() or x == "."])
        v = float(v)
        if v < cutoff:
            return None
        else:
            return 0
    if operator == "<":
        v = "".join([x for x in str(v) if x.isdigit() or x == "."])
        v = float(v)
        if v > cutoff:
            return None
        else:
            return 1

mic_cols=[]
for c in CUTOFFS_CONC:
    all_mic[f"mic_{c}um"] = all_mic["mic_um"].apply(lambda x: binarizer(x, c))
    mic_cols += [f"mic_{c}um"]

#STEP 7: Merge duplicates using Binary data

def aggregate_and_binarize(df, smiles_col, mic_cols, threshold=0.5):
    aggregated = df.groupby(smiles_col)[mic_cols].mean()
    for col in mic_cols:
        aggregated[col] = (aggregated[col] > threshold).astype(int)
    aggregated = aggregated.reset_index()
    return aggregated

all_mic = aggregate_and_binarize(all_mic, "smiles", mic_cols)
print(all_mic.shape)
all_mic.to_csv(os.path.join(processed_dir, f"{pathogen_code}_spark_mic.csv"), index=False)

"""
df = pd.read_csv(os.path.join(data_dir, 'spark', 'SPARK MIC Data.csv'))
print(df.shape)

# STEP 1: Eliminate rows without SMILES
df = df[df["SMILES"].notnull()]
print(df.shape)

# STEP 2: Separate Curated vs Extracted
df_curated = df[['SMILES'] + df.columns[df.columns.str.contains("Curated & Transformed")].tolist()]
print(df_curated.shape)
df_extracted = df[['SMILES'] + df.columns[df.columns.str.contains("Extracted & Uploaded")].tolist()]
print(df_extracted.shape)

# STEP 3: in each dataset, keep the pathogen of interest
df_curated = df_curated[df_curated["Curated & Transformed MIC Data: Species"]==pathogen_name]
print(df_curated.shape)
df_extracted = df_extracted[df_extracted["Extracted & Uploaded MIC Data: Species"]==pathogen_name]
print(df_extracted.shape)

# STEP 4: in each dataset, consolidate the experiments of interest
df_curated_mic = df_curated[~df_curated["Curated & Transformed MIC Data: MIC (in µM) (µM)"].isna()]
df_curated_gic = df_curated[~df_curated["Curated & Transformed MIC Data: GIC % inhibition"].isna()]
print(len(df_curated_mic), len(df_curated_gic))

df_extracted_mic = df_extracted[~df_extracted["Extracted & Uploaded MIC Data: MIC (in µM) (µM)"].isna()]
df_extracted_gic = df_extracted[~df_extracted["Extracted & Uploaded MIC Data: GIC % Inhibition (%)"].isna()]
print(len(df_extracted_mic), len(df_extracted_gic))

# STEP 5: Eliminate SMILES which cannot be converted
df_curated_mic = standardise_smiles(df_curated_mic)
df_curated_gic = standardise_smiles(df_curated_gic)
df_extracted_mic = standardise_smiles(df_extracted_mic)
df_extracted_gic = standardise_smiles(df_extracted_gic)

print(len(df_curated_mic), len(df_curated_gic), len(df_extracted_mic), len(df_extracted_gic))

# STEP 6: convert units to binary
df_curated_mic = df_curated_mic[["std_smiles", "Curated & Transformed MIC Data: MIC (in µM) (µM)"]]
df_extracted_mic = df_extracted_mic[["std_smiles", "Extracted & Uploaded MIC Data: MIC (in µM) (µM)"]]
df_curated_mic.rename(columns={"Curated & Transformed MIC Data: MIC (in µM) (µM)":"mic_um", "std_smiles": "smiles"}, inplace=True)
df_extracted_mic.rename(columns={"Extracted & Uploaded MIC Data: MIC (in µM) (µM)":"mic_um", "std_smiles": "smiles"}, inplace=True)

df_curated_mic.to_csv(os.path.join(data_dir, "spark", "curated_mic_wip.csv"), index=False)

def extract_operator(v):
    v = str(v)
    try:
        float(v)
        return "="
    except:
        if v.startswith(">"):
            return ">"
        elif v.startswith("<"):
            return "<"
        return None


def binarizer(v, cutoff):
    operator = extract_operator(v)
    if operator is None:
        return None
    if operator == "=":
        v = float(v)
        if v <= cutoff:
            return 1
        else:
            return 0
    if operator == ">":
        v = "".join([x for x in str(v) if x.isdigit() or x == "."])
        v = float(v)
        if v < cutoff:
            return None
        else:
            return 0
    if operator == "<":
        v = "".join([x for x in str(v) if x.isdigit() or x == "."])
        v = float(v)
        if v > cutoff:
            return None
        else:
            return 1

mic_cols=[]
for c in CUTOFFS_CONC:
    df_curated_mic[f"mic_{c}um"] = df_curated_mic["mic_um"].apply(lambda x: binarizer(x, c))
    df_extracted_mic[f"mic_{c}um"] = df_extracted_mic["mic_um"].apply(lambda x: binarizer(x, c))
    mic_cols += [f"mic_{c}um"]

print(df_curated_mic.shape, df_extracted_mic.shape)

# STEP 7: merge duplicates
print(df_curated_mic["smiles"].nunique(), df_extracted_mic["smiles"].nunique())

df_mic = pd.concat([df_curated_mic, df_extracted_mic], ignore_index=True)

def aggregate_and_binarize(df, smiles_col, mic_cols, threshold=0.5):
    # Group by the smiles column and calculate the mean for mic columns
    aggregated = df.groupby(smiles_col)[mic_cols].mean()
    # Convert averages to binary based on the threshold
    for col in mic_cols:
        aggregated[col] = (aggregated[col] > threshold).astype(int)
    # Reset index to include smiles as a column
    aggregated = aggregated.reset_index()
    return aggregated

df_mic = aggregate_and_binarize(df_mic, "smiles", mic_cols)
print(df_mic.shape)
df_mic.to_csv(os.path.join(processed_dir, f"{pathogen_code}_spark_mic.csv"), index=False)
"""
