import glob
from rdkit.Chem import PandasTools
import pandas as pd


scored_sdfs = glob.glob('../files_to_analyse/*_scored.sdf')
dfs = []
for f in scored_sdfs:
    df = PandasTools.LoadSDF(f)
    dfs.append(df)

all_mols = pd.concat(dfs)
props = [key for key in all_mols.keys() if key not in ['ID', 'ROMol']]
PandasTools.WriteSDF(all_mols, '../final_scored_agg.sdf', properties=props)
all_mols.to_csv('../final_scored_agg.csv', columns=props)