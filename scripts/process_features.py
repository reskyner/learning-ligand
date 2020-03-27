"""Cleans and processes pre-computed RDKit, RF-Score, and BINANA features.

Usage:
    process_features.py <pdbbind_data_file> [<core_set_file>]

Arguments:
    pdbbind_index_file INDEX_general_PL_data.2018 (included with PDBbind, not redistributed here)  
    core_set_file      optional file containing PDB ids of core set structures

Loads the features computed by:
    compute_rdkit_features.py,
    compute_rfscore_features.py,
    compute_binana_features.py,
and the binding data provided with the PDBbind 2018 general set. Samples
for which features cannot be computed are dropped from all feature sets.

Drops features with zero variance, ignoring any PDB codes provided in the
(optional) PDBbind core set (or any other list of PDB codes).

Each processed feature set is saved in .csv format, and if a version of the
core set is provided, the PDB codes corresponding to the training and test set
used when removing zero-variance features are also saved to .txt files.

"""

import pathlib

import numpy as np
import pandas as pd

from docopt import docopt

# process command line arguments
args = docopt(__doc__)
pdbbind_data_file = args['<pdbbind_data_file>']
core_set_file = args['<core_set_file>']

# load computed features
rdkit_file = pathlib.Path('..', 'data', 'pdbbind_2018_general_rdkit_features.csv')
rfscore_file = pathlib.Path('..', 'data', 'pdbbind_2018_general_rfscore_features.csv')
nnscore_file = pathlib.Path('..', 'data', 'pdbbind_2018_general_binana_features.csv')

rdkit_features = pd.read_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_rdkit_features.csv'), index_col=0)
rfscore_features = pd.read_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_rfscore_features.csv'), index_col=0)
nnscore_features = pd.read_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_binana_features.csv'), index_col=0)

# Drop 'Ipc' from RDKit feature set
# It generates very large values for larger molecules. https://github.com/rdkit/rdkit/issues/1527
rdkit_features = rdkit_features.drop(['Ipc'], axis='columns')

# load PDBbind binding data
with open(pdbbind_data_file) as f:
    pdbbind_data = [line.strip().split() for line in f if not line.startswith('#')]
binding_data = {line[0]: float(line[3]) for line in pdbbind_data}
binding_data = pd.DataFrame.from_dict(binding_data, orient='index')

# Drop any data points will null-valued descriptors
rdkit_features = rdkit_features.replace([np.inf, -np.inf], np.nan)
rdkit_features = rdkit_features.dropna(axis='index', how='any')
rfscore_features = rfscore_features.dropna(axis='index', how='any')
nnscore_features = nnscore_features.dropna(axis='index', how='any')

# Keep only data for which allfeatures were computed
pdbs_to_use = nnscore_features.index.intersection(rdkit_features.index.intersection(rfscore_features.index))
binding_data = binding_data.loc[pdbs_to_use]
rdkit_features = rdkit_features.loc[pdbs_to_use]
rfscore_features = rfscore_features.loc[pdbs_to_use]
nnscore_features = nnscore_features.loc[pdbs_to_use]

rfscore_features = rfscore_features.sort_index()
nnscore_Features = nnscore_features.sort_index()
rdkit_features = rdkit_features.sort_index()
binding_data = binding_data.sort_index()

# if specified, exclude core set data when removing zero-variance features
if core_set_file:
    with open(core_set_file, 'r') as f:
        core_set = [line.strip() for line in f]
    test_set = [pdb for pdb in core_set if pdb in pdbs_to_use]
    training_set = [pdb for pdb in pdbs_to_use if pdb not in test_set]
else:
    training_set = pdbs_to_use

# Drop features with zero variance across the training set
to_drop = rdkit_features.columns[(rdkit_features.loc[training_set].var() == 0)]
rdkit features = rdkit_features.drop(to_drop, axis='columns')
to_drop = rfscore_features.columns[(rfscore_features.loc[training_set].var() == 0)]
rfscore_features = rfscore_features.drop(to_drop, axis='columns')
to_drop = nnscore_features.columns[(nnscore_features.loc[training_set].var() == 0)]
nnscore_Features = nnscore_features.drop(to_drop, axis='columns')

# if a core set was excluded from the training data, save a list of
# PDB codes to use as trainining and test sets
if core_set_file:
    with open(os.path.join('..', 'data', 'training_set.txt'), 'w') as f:
        print(*training_set, sep='\n', file=f)
    with open(os.path.join('..', 'data', 'test_set.txt'), 'w') as f:
        print(*test_set, sep='\n', file=f)

# save processed feature sets and binding data
rdkit_features.to_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_rdkit_features_clean.csv'))
rfscore_features.to_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_rfscore_features_clean.csv'))
nnscore_features.to_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_binana_features_clean.csv'))
binding_data.to_csv(pathlib.Path('..', 'data', 'pdbbind_2018_general_binding_data_clean.csv'))

