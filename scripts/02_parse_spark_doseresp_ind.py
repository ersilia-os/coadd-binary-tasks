import os
import pandas as pd
import sys

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from utils import binarizer, merge_replicas
from spark_config import spark_data

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, '..', 'data'))
processed_dir = os.path.abspath(os.path.join(root, '..', 'processed'))
config_dir = os.path.abspath(os.path.join(root, '..', 'config'))

pathogen_code = "abaumannii"
dc = pd.read_csv(os.path.join(config_dir, 'pathogens.csv'))
pathogen_name = dc[dc['pathogen_code'] == pathogen_code]['search_text'].values[0]

CUTOFFS_CONC = [1,10,40]
cutoff_cols = [f"mic_{c}um" for c in CUTOFFS_CONC]

print("Working on pathogen: ", pathogen_name, pathogen_code)

print("Starting with primary screening data")

spark_data = spark_data[pathogen_code]

for k,v in spark_data.items():
    print("--------------------")
    print(f"Working with dataset contributed by {k}")
    df = pd.read_csv(os.path.join(data_dir, "spark", "preprocessed", f"spark_{k}.csv"), low_memory=False)
    values_col = v["value_col"]
    species_col = v["sp_col"]
    strain_col = v["strain_col"]
    strain = v["strain"]
    df = df[df[species_col]==pathogen_name]
    if strain != None:
        df = df[df[strain_col]==strain]
    df = df[["smiles", values_col]]
    for c in CUTOFFS_CONC:
            df[f"mic_{c}um"] = df[values_col].apply(lambda x: binarizer(x, c))
    rows_with_none = df[df.isnull().any(axis=1)]
    print("SMILES that cannot be binarised at one or more cut-offs: ", len(rows_with_none))
    df_dicts = []
    for c in cutoff_cols:
        result = merge_replicas(df, c)
        df_dicts += [result]

    smiles_lists = [d["smiles"] for d in df_dicts]
    if not all(s == smiles_lists[0] for s in smiles_lists):
        raise ValueError("The SMILES differ between cut-offs, please check")
    df = pd.DataFrame({"smiles": smiles_lists[0]})

    for d in df_dicts:
        for key, values in d.items():
            if key != "smiles":
                df[key] = values
    print("Total datapoints: ", len(df))
    df.to_csv(os.path.join(processed_dir, "spark_individual", f"abaumannii_spark_doseresponse_{k}.csv"), index=False)


#Join all files
dfs = []
for k,v in spark_data.items():
    df = pd.read_csv(os.path.join(processed_dir, "spark_individual", f"{pathogen_code}_spark_doseresponse_{k}.csv"))
    dfs += [df]
df = pd.concat(dfs)
print("--------------------")
print("Merging all individual files")
print("All datapoints in individual files merged: ", len(df))
print("Unique datapoints in individual files merged: ", len(set(df["smiles"])))
df_dicts = []
for c in cutoff_cols:
    result = merge_replicas(df, c)
    df_dicts += [result]

smiles_lists = [d["smiles"] for d in df_dicts]
if not all(s == smiles_lists[0] for s in smiles_lists):
    raise ValueError("The SMILES differ between cut-offs, please check")
df = pd.DataFrame({"smiles": smiles_lists[0]})

for d in df_dicts:
    for key, values in d.items():
        if key != "smiles":
            df[key] = values
print("Total datapoints: ", len(df))
df.to_csv(os.path.join(processed_dir, "spark_individual", f"{pathogen_code}_spark_doseresponse_merged.csv"), index=False)