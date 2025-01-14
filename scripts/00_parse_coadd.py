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

CUTOFFS = [50, 80, 90]

print("Working on pathogen: ", pathogen_name, pathogen_code)

print("Starting with primary screening data")

# Read in the coadd data

df = pd.read_csv(os.path.join(data_dir, 'coadd', 'CO-ADD_InhibitionData_r03_01-02-2020_CSV.csv'))
df = df[df["INHIB_AVE"].notnull()]
df = df[df["SMILES"].notnull()]
df = df[df["INHIB_AVE"] > -50]
print(df["CONC"].value_counts(ascending=False))

df = df[df['ORGANISM'] == pathogen_name]

# Standardise the SMILES
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

data = {
    "low_conc": df[df['CONC'] != '300 ug/mL'][['std_smiles', 'INHIB_AVE']].rename(columns={"INHIB_AVE": "inhib_perc"}),
    "high_conc": df[df['CONC'] == '300 ug/mL'][['std_smiles', 'INHIB_AVE']].rename(columns={"INHIB_AVE": "inhib_perc"})
}

binary_data = collections.OrderedDict()
for k, v in data.items():
    data_ = collections.OrderedDict()
    data_["smiles"] = list(v['std_smiles'])
    for cutoff in CUTOFFS:
        inhib_list = list(v['inhib_perc'])
        inhib_list = [1 if x >= cutoff else 0 for x in inhib_list]
        data_["inhib_{0}".format(cutoff)] = inhib_list
    data_ = pd.DataFrame(data_)
    data_ = data_[data_["smiles"].notnull()]
    data_.to_csv(os.path.join(processed_dir, '{0}_coadd_inhib_{1}.csv'.format(pathogen_code, k)), index=False)


print("Finished with primary screening data")
print("Starting with secondary screening data")

df = pd.read_csv(os.path.join(data_dir, 'coadd', 'CO-ADD_DoseResponseData_r03_01-02-2020_CSV.csv'))
print(df[df["DRVAL_TYPE"].isnull()])
print(df["DRVAL_TYPE"].value_counts())

df = df[df['ORGANISM'] == pathogen_name]
df = df[df['DRVAL_TYPE'] == 'MIC']

print("Number of compounds: ", df.shape[0])

def convert_ugml_to_uM(smiles, concentration_ugml):
    """
    Convert a concentration from µg/mL to µM using the molecular weight from a SMILES string.

    Parameters:
    - smiles (str): SMILES string of the compound.
    - concentration_ugml (float): Concentration in µg/mL.

    Returns:
    - float: Concentration in µM.
    """
    try:
        # Create a molecule object from the SMILES string
        molecule = Chem.MolFromSmiles(smiles)
        if not molecule:
            raise ValueError("Invalid SMILES string")
        
        # Calculate molecular weight
        molecular_weight = Descriptors.MolWt(molecule)
        
        # µg/mL to µM conversion formula
        # µM = (concentration in µg/mL) / (molecular weight in g/mol) * 1000
        concentration_uM = (concentration_ugml / molecular_weight) * 1000
        return concentration_uM
    except Exception as e:
        print(f"Error: {e}")
        return None

df = df[df["DRVAL_UNIT"].notnull()]

R = []
for v in df[["SMILES", "DRVAL_UNIT", "DRVAL_MEDIAN", "DRVAL_TYPE"]].values:
    smiles = v[0]
    unit = v[1]
    value = str(v[2])
    valtype = v[3]
    if value.startswith("<="):
        value = value[2:]
        operator = "<"
    elif value.startswith("<"):
        value = value[1:]
        operator = "<"
    elif value.startswith(">="):
        value = value[2:]
        operator = ">"
    elif value.startswith(">"):
        value = value[1:]
        operator = ">"
    elif value.startswith("="):
        value = value[1:]
        operator = "="
    else:
        operator = "="
    value = float(value)
    if unit == "uM":
        R.append([value, operator, valtype])
    elif unit == "ug/mL":
        R.append([convert_ugml_to_uM(smiles, value), operator, valtype])
    else:
        raise Exception("Unknown unit", unit)
    
dp = pd.DataFrame(R, columns=["value_um", "operator", "type"])
std_smiles = []
for smiles in tqdm(list(df["SMILES"])):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = standardise.run(mol)
        std_smiles.append(Chem.MolToSmiles(mol))
    except:
        std_smiles.append(None)

dp["smiles"] = std_smiles
dp = dp[dp["smiles"].notnull()]

CUTOFFS = [50, 100]

def get_binary_dataset(dp, cutoff):
    bin_values = []
    smiles_list = []
    for v in dp.values:
        value = v[0]
        operator = v[1]
        smiles = v[3]
        if operator == ">" and value < cutoff:
            continue
        if operator == "<" and value > cutoff:
            continue
        if operator == ">":
            bin_values.append(0)
            smiles_list.append(smiles)
            continue
        if operator == "<":
            bin_values.append(1)
            smiles_list.append(smiles)
            continue
        if operator == "=":
            if value <= cutoff:
                bin_values.append(1)
                smiles_list.append(smiles)
            else:
                bin_values.append(0)
                smiles_list.append(smiles)
            continue
    dr = pd.DataFrame({"smiles": smiles_list, "mic_{0}um".format(cutoff): bin_values})
    return dr

for cutoff in CUTOFFS:
    dr = get_binary_dataset(dp, cutoff)
    dr.to_csv(os.path.join(processed_dir, pathogen_code + f"_doseresp_{cutoff}um.csv"), index=False)