import os
import pandas as pd
import collections

import sys

abspath = os.path.abspath(__file__)
sys.path.append(abspath)

from utils import standardise_smiles
from spark_config import spark_data

root = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(os.path.join(root, '..', 'data'))
processed_dir = os.path.abspath(os.path.join(root, '..', 'processed'))
config_dir = os.path.abspath(os.path.join(root, '..', 'config'))

pathogen_code = "abaumannii"
dc = pd.read_csv(os.path.join(config_dir, 'pathogens.csv'))
pathogen_name = dc[dc['pathogen_code'] == pathogen_code]['search_text'].values[0]


df = pd.read_csv(os.path.join(data_dir, "spark", "preprocessed", f"spark_co-add.csv"), low_memory=False)
val_col = 'CO-ADD Single Concentration Data (Inhibition): Inhibition % (%)'
conc_col = 'CO-ADD Single Concentration Data (Inhibition): Inhibition Concentration (µg/mL)'
species_col = 'CO-ADD Single Concentration Data (Inhibition): Species'

CUTOFFS = [50, 80, 90]

df = df[df[species_col]==pathogen_name]
df = df[df[val_col].notna()]

print(f"Total Datapoints with Single point inhibition data in {pathogen_code}: ", len(set(df["smiles"])))
print(f"Total unique SMILES with Single point inhibition data in {pathogen_code}: ", len(set(df["smiles"])))

data = {
    "low_conc": df[df[conc_col] < 300][['smiles', val_col]].rename(columns={val_col: "inhib_perc"}),
    "high_conc": df[df[conc_col] == 300 ][['smiles', val_col]].rename(columns={val_col: "inhib_perc"})
}


binary_data = collections.OrderedDict()
for k, v in data.items():
    data_ = collections.OrderedDict()
    data_["smiles"] = list(v['smiles'])
    for cutoff in CUTOFFS:
        inhib_list = list(v["inhib_perc"])
        inhib_list = [1 if x >= cutoff else 0 for x in inhib_list]
        data_["inhib_{0}".format(cutoff)] = inhib_list
    data_ = pd.DataFrame(data_)
    data_ = data_[data_["smiles"].notnull()]
    data_.to_csv(os.path.join(processed_dir, f'{pathogen_code}_spark_coadd_inhib_{k}.csv'), index=False)
    print(f"{k} total data: ", len(set(data_["smiles"])), "Positives 50: ", len(data_[data_["inhib_50"] == 1]), "Positives 80: ", len(data_[data_["inhib_80"] == 1]), "Positives 90: ", len(data_[data_["inhib_90"] == 1]))



coadd = pd.read_csv(os.path.join(processed_dir, 'abaumannii_coadd_inhib_low_conc.csv'))










"""

df = standardise_smiles(df, "smiles")
coadd = standardise_smiles(coadd, "smiles")

print("Difference between Spark and CoAdd: ", len(set(df["smiles"]) - set(coadd["smiles"])))
print("Difference between CoaDD and Spark: ", len(set(coadd["smiles"]) - set(df["smiles"])))

data = {
    "low_conc": df[df['CONC'] != '300 ug/mL'][['std_smiles', 'INHIB_AVE']].rename(columns={"INHIB_AVE": "inhib_perc"}),
    "high_conc": df[df['CONC'] == '300 ug/mL'][['std_smiles', 'INHIB_AVE']].rename(columns={"INHIB_AVE": "inhib_perc"})
}
"""