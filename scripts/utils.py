from tqdm import tqdm
from rdkit import Chem
import numpy as np
from standardiser import standardise
import collections

def standardise_smiles(df, smi_col):
    std_smiles = []
    for smiles in tqdm(list(df[smi_col])):
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = standardise.run(mol)
            std_smiles.append(Chem.MolToSmiles(mol))
        except:
            std_smiles.append(None)
    df["std_smiles"] = std_smiles
    df = df[df["std_smiles"].notnull()]
    return df

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

def aggregate_bin_cols(df, smiles_col, mic_cols, threshold=0.5):
    aggregated = df.groupby(smiles_col)[mic_cols].mean()
    for col in mic_cols:
        aggregated[col] = (aggregated[col] >= threshold).astype(int)
    aggregated = aggregated.reset_index()
    return aggregated


def merge_replicas(data, col, threshold=0.5): #TODO add a list of SMILES with over 5 replicas and delete them?
    merged = collections.defaultdict(list)
    smi_list=data["smiles"].tolist()
    val_list = data[col].tolist()
    for i,s in enumerate(smi_list):
        val = val_list[i]
        merged[s].append(val)
    result = collections.defaultdict(list)
    for s, val in merged.items():
        valid_values = [v for v in val if v is not None and not np.isnan(v)] 
        if valid_values:
            mean_value = sum(valid_values) / len(valid_values)
            final_value = 1 if mean_value >= threshold else 0
        else:
            final_value = None # If no valid values, keep None
        result["smiles"].append(s)
        result[col].append(final_value)
    return result