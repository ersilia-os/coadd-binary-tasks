import pandas as pd
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)

from ..src.utils import standardise_smiles

data_dir = os.path.abspath(os.path.join(root, '..', 'data', 'spark'))

spark_files = {'SPARK Accumulation Data.csv': 'spark_accumulation.csv', 
               'SPARK Data Merck & Kyorin Contribution.csv': 'spark_merck.csv', 
               'SPARK Data Quave Lab {Emory University} Publications.csv': 'spark_quave.csv', 
               'SPARK Data Compounds & Physicochemical Properties.csv': 'spark_physicochemical.csv', 
               'SPARK MIC Data.csv': 'spark_mic.csv', 
               'SPARK IC50 Data.csv': 'spark_ic50.csv', 
               'SPARK Data CO-ADD Contribution.csv': 'spark_co-add.csv',
               'SPARK Data Achaogen Contribution.csv': 'spark_achaogen.csv', 
               'SPARK Data Novartis Contribution.csv': 'spark_novartis.csv'}


for k,v in spark_files.items():
    print(k)
    df = pd.read_csv(os.path.join(data_dir, k))
    df1 = df[df["SMILES"].notna()]
    df2 = standardise_smiles(df1, "SMILES")
    df2 = df2.drop(columns=["SMILES", " SPARK Data Downloads"]+[col for col in df2.columns if "Terms of Use" in col])
    df2 = df2.rename(columns={"std_smiles": "smiles"})
    df2 = df2.dropna(axis=1, how='all')
    print("Initial compounds: ", len(df), "No SMILES: ", len(df) - len(df1), "Non standardised: ", len(df1)-len(df2))
    df2.to_csv(os.path.join(data_dir, "preprocessed",v), index=False)