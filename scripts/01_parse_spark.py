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

UTOFFS = [50, 80, 90]

print("Working on pathogen: ", pathogen_name, pathogen_code)

print("Starting with primary screening data")

# Read in the spark data (MIC)

df = pd.read_csv(os.path.join(data_dir, 'spark', 'SPARK MIC Data.csv'))
df = df[df["SMILES"].notnull()]
df = df[df['ORGANISM'] == pathogen_name]