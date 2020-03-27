from rdkit.Chem import PandasTools
from compute_rfscore_features import featurise_multi_ligs
from compute_rdkit_features import rdkit_features_multilig
from compute_binana_features import binana_multi_ligs
import glob
import pandas as pd

lines = []

# clean aggregated file
for line in open('/Users/uzw12877/Code/ligand-learn/files_to_analyse/Mpro-17-aggregated-best-transfs.sdf').readlines():
    if '$>' in line:
        new_line = line.replace('$>','>')
        print(new_line)
        lines.append(new_line)
    else:
        lines.append(line)

with open('files_to_analyse/cleaned.sdf', 'a') as f:
    for line in lines:
        f.write(line)
f.close()

# split agreggated file into one sdf for each protein
df = PandasTools.LoadSDF('../files_to_analyse/cleaned.sdf')
props = [key for key in df.keys() if key not in ['ID', 'ROMol']]
targets = list(set(df['Protein_Target']))
print(targets)
sdfs_dict = {}
for target in targets:
    sdfs_dict[target] = df[df['Protein_Target']==target]

for key in sdfs_dict.keys():
    PandasTools.WriteSDF(sdfs_dict[key], f'../files_to_analyse/{key}_mols.sdf', properties=props)

# calculate rf-score features
for f in glob.glob('../files_to_analyse/Mpro-x*.sdf'):
    name = f.split('/')[-1].replace('.sdf', '').replace('_mols', '')
    protein = f'../files_to_analyse/prots/{name}_apo-desolv.pdb'
    print(protein)

    results = featurise_multi_ligs(protein, f)

    features = pd.DataFrame.from_dict(results).T
    print(features)
    features.to_csv(f'../files_to_analyse/features/rfscore-features_{name}.csv')

# calculate rdkit features
for f in glob.glob('../files_to_analyse/Mpro-x*.sdf'):
    name = f.split('/')[-1].replace('.sdf', '').replace('_mols', '')
    features = rdkit_features_multilig(f)

    features = pd.DataFrame.from_dict(features).T
    features.to_csv(f'../files_to_analyse/features/rdkit-features_{name}.csv')

# calculate binana features
for f in glob.glob('../files_to_analyse/Mpro-x*.sdf'):
    name = f.split('/')[-1].replace('.sdf', '').replace('_mols', '')
    protein = f'../files_to_analyse/prots/{name}_apo-desolv.pdb'
    print(protein)

    results = binana_multi_ligs(protein, f)

    features = pd.DataFrame.from_dict(results).T
    print(features)
    features.to_csv(f'../files_to_analyse/features/binana-features_{name}.csv')
