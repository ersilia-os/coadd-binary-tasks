import os
import pandas as pd

root = os.path.join(os.path.dirname(__file__))
data_dir = os.path.abspath(os.path.join(root, '..', 'data'))
processed_dir = os.path.abspath(os.path.join(root, '..', 'processed'))

# Read in the coadd data

df = pd.read_csv(os.path.join(data_dir, 'CO-ADD_InhibitionData_r03_01-02-2020_CSV.csv'))

print(df.head())