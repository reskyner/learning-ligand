"""
Quick and dirty clustering of ligands in the PDBbind 2018 general set 
at Morgan FP Tanimoto similarity > 0.9 with respect to ligands in the 
2007, 2013, and 2016 core sets.
"""

import pandas as pd

# previously computed TC between core set and general set ligands
tc_file = '../data/pdbbind_2018_general_ligand_tanimoto_coefficients.csv'
tc = pd.read_csv(tc_file, index_col=0)

with open('../data/pdbbind_2007_core_pdbs.txt') as f:
    core_2007 = [l.strip() for l in f if l.strip() in tc.index]

with open('../data/pdbbind_2013_core_pdbs.txt') as f:
    core_2013 = [l.strip() for l in f if l.strip() in tc.index]

with open('../data/pdbbind_2016_core_pdbs.txt') as f:
    core_2016 = [l.strip() for l in f if l.strip() in tc.index]

core_2007_similar = []

for pdb in core_2007:
    core_2007_similar.extend([pdb for pdb in tc[tc[pdb]>=0.9].index])
core_2007_similar = sorted(list(set(core_2007_similar)))

core_2013_similar = []

for pdb in core_2013:
    core_2013_similar.extend([pdb for pdb in tc[tc[pdb]>=0.9].index])
core_2013_similar = sorted(list(set(core_2013_similar)))

core_2016_similar = []

for pdb in core_2016:
    core_2016_similar.extend([pdb for pdb in tc[tc[pdb]>=0.9].index])
core_2016_similar = sorted(list(set(core_2016_similar)))

print(len(core_2007_similar), len(core_2013_similar), len(core_2016_similar))

with open('../data/pdbbind_2018_general_2007_core_similar_ligands.txt', 'w') as f:
    print(*core_2007_similar, sep='\n', file=f)
with open('../data/pdbbind_2018_general_2013_core_similar_ligands.txt', 'w') as f:
    print(*core_2013_similar, sep='\n', file=f)
with open('../data/pdbbind_2018_general_2016_core_similar_ligands.txt', 'w') as f:
    print(*core_2016_similar, sep='\n', file=f)

