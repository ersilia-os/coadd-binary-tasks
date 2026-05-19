from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
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
    """Return the comparison operator prefix of a censored value string.

    Handles >=, <=, >, <, = prefixes and bare numbers (treated as =).
    Returns one of '=', '>', '<', or None if the string is unparseable.
    Also returns the remainder string with the operator stripped.

    Examples:
        '>64'   -> ('>', '64')
        '<=0.5' -> ('<', '0.5')
        '32'    -> ('=', '32')
    """
    s = str(v).strip()
    if s.startswith(">="):
        return ">", s[2:]
    if s.startswith("<="):
        return "<", s[2:]
    if s.startswith(">"):
        return ">", s[1:]
    if s.startswith("<"):
        return "<", s[1:]
    if s.startswith("="):
        return "=", s[1:]
    return "=", s


def parse_numeric(val_str):
    """Parse a censored numeric string and return (operator, value).

    Returns ("=", float) for plain numbers, (">", float) / ("<", float) for
    censored values, and (None, None) if the numeric part cannot be parsed.

    Examples:
        '>64'   -> ('>', 64.0)
        '<=0.5' -> ('<', 0.5)
        '32'    -> ('=', 32.0)
    """
    op, s = extract_operator(val_str)
    if op is None:
        return None, None
    try:
        return op, float(s)
    except ValueError:
        return None, None

def binarizer(v, cutoff):
    operator, s = extract_operator(v)
    if operator is None:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if operator == "=":
        return 1 if v <= cutoff else 0
    if operator == ">":
        return None if v < cutoff else 0
    if operator == "<":
        return None if v > cutoff else 1

def aggregate_bin_cols(df, smiles_col, mic_cols, threshold=0.5):
    aggregated = df.groupby(smiles_col)[mic_cols].mean()
    for col in mic_cols:
        aggregated[col] = (aggregated[col] >= threshold).astype(int)
    aggregated = aggregated.reset_index()
    return aggregated


def convert_ugml_to_uM(smiles, concentration_ugml):
    """Convert concentration from µg/mL to µM using molecular weight from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        mw = Descriptors.MolWt(mol)
        return (concentration_ugml / mw) * 1000
    except Exception:
        return None


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