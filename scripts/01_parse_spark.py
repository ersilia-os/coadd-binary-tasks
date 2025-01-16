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
    for smiles in tqdm(list(df["SMILES"])):
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = standardise.run(mol)
            std_smiles.append(Chem.MolToSmiles(mol))
        except:
            std_smiles.append(None)
    print(len(std_smiles), df.shape)
    df["std_smiles"] = std_smiles
    df = df[df["std_smiles"].notnull()]
    return df

CUTOFFS_PERC = [50, 80, 90]
CUTOFFS_CONC = [10,25,50]

print("Working on pathogen: ", pathogen_name, pathogen_code)

print("Starting with primary screening data")

# Read in the spark data (MIC)
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
cols = df_curated_mic.columns

# STEP 5: Eliminate SMILES which cannot be converted
df_curated_mic = standardise_smiles(df_curated_mic)
df_curated_gic = standardise_smiles(df_curated_gic)
df_extracted_mic = standardise_smiles(df_extracted_mic)
df_extracted_gic = standardise_smiles(df_extracted_gic)

print(len(df_curated_mic), len(df_curated_gic), len(df_extracted_mic), len(df_extracted_gic))

# STEP 6: convert units to binary
df_curated_mic = df_curated_mic[["std_smiles", "Curated & Transformed MIC Data: MIC (in µM) (µM)"]]
df_extracted_mic = df_extracted_mic[["std_smiles", "Extracted & Uploaded MIC Data: MIC (in µM) (µM)"]]
df_curated_mic.rename(columns={"Curated & Transformed MIC Data: MIC (in µM) (µM)":"MIC", "std_smiles": "smiles"}, inplace=True)
df_extracted_mic.rename(columns={"Extracted & Uploaded MIC Data: MIC (in µM) (µM)":"MIC", "std_smiles": "smiles"}, inplace=True)

def convert_to_binary(df, cutoffs, colname):
    valid_values=[]
    for item in df[colname]:
        try:
           value = float(item)
        except:
            print(item)
        valid_values += [value]
    for cutoff in cutoffs:
        bin = [1 if x <= cutoff else 0 for x in valid_values]
        df[f"{colname}_binary_{cutoff}"] = bin
    return df

df_curated_mic = convert_to_binary(df_curated_mic, CUTOFFS_CONC, "MIC")
df_extracted_mic = convert_to_binary(df_extracted_mic, CUTOFFS_CONC, "MIC")

print(df_curated_mic.shape, df_extracted_mic.shape)

# STEP 7: merge duplicates
