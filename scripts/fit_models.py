"""Fits RF models and saves affinity predictions for PDBbind core sets.

"""
import itertools
import json
import os
import joblib
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestRegressor
from scipy.stats import pearsonr

print('Loading data')

rdkit_features = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_rdkit_features_clean.csv'), index_col=0)
rfscore_features = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_rfscore_features_clean.csv'), index_col=0)
nnscore_features = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_binana_features_clean.csv'), index_col=0)
binding_data = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_binding_data_clean.csv'), index_col=0, squeeze=True)
binding_data = binding_data.rename('pK')

all_features = pd.concat([rdkit_features, rfscore_features, nnscore_features], axis='columns')

feature_sets = {
    'RDKit': rdkit_features,
    'RF-Score': rfscore_features,
    'RF-Score + RDKit': pd.concat([rdkit_features, rfscore_features], axis='columns'),
    'NNScore 2.0': nnscore_features,
    'NNScore 2.0 + RDKit': pd.concat([rdkit_features, nnscore_features], axis='columns')

}

vina_terms = pd.Index(['vina_gauss1', 'vina_gauss2', 'vina_hydrogen', 'vina_hydrophobic', 'vina_repulsion', 'num_rotors'])
feature_sets['Vina'] = all_features[vina_terms]
feature_sets['Vina + RDKit'] = all_features[vina_terms.union(rdkit_features.columns)]
feature_sets['Vina + RDKit'] = feature_sets['Vina + RDKit'].drop('NumRotatableBonds', axis='columns')
feature_sets['RF-Score v3'] = all_features[vina_terms.union(rfscore_features.columns)]
feature_sets['RF-Score v3 + RDKit'] = pd.concat([feature_sets['RF-Score v3'], feature_sets['RDKit']], axis='columns')

print('Clustering data')

# Load general and refined sets
data_sets = {}
for year, subset in itertools.product([2013, 2014, 2015, 2016, 2017, 2018], ['general', 'refined']):
    with open(f'../data/pdbbind_{year}_{subset}_pdbs.txt') as f:
        data_sets[f'{year} {subset}'] = [l.strip() for l in f]

# Load core sets
core_sets = {}
for year in ['2007', '2013', '2016']:
    with open(f'../data/pdbbind_{year}_core_pdbs.txt') as f:
        core_sets[year] = sorted([l.strip() for l in f])
core_sets['all'] = [pdb for pdb in core_sets['2007']]
core_sets['all'] = core_sets['all'] + [pdb for pdb in core_sets['2013'] if pdb not in core_sets['all']]
core_sets['all'] = core_sets['all'] + [pdb for pdb in core_sets['2016'] if pdb not in core_sets['all']]


# load BLASTclust clustering of the PDB
blast_clusters = {}
for cutoff in [30, 40, 50, 70, 90, 95, 100]:
    with open(f'../data/bc-{cutoff}.out') as f:
        blast_clusters[cutoff] = [set(item[:4].lower() for item in line.strip().split()) for line in f]

# list PDBbind structures in the same BLAST cluster as any of the core set structures
core_blast_similar_pdbs = {}
for core_set in core_sets:
    core_blast_similar_pdbs[core_set] = {}
    for cutoff in blast_clusters:
        core_blast_similar_pdbs[core_set][cutoff] = set()
        for pdb in core_sets[core_set]:
            for cluster in blast_clusters[cutoff]:
                if pdb in cluster:
                    core_blast_similar_pdbs[core_set][cutoff].update(cluster)

# for each release of PDBbind, list the structures in the general and refined set that belong to the same BLAST cluster as a core set structure
overlap_size_general = {}
overlap_size_refined = {}
overlap_general = {}
overlap_refined = {}

for core_set in core_sets:
    overlap_size_general[core_set] = {}
    overlap_size_refined[core_set] = {}
    overlap_general[core_set] = {}
    overlap_refined[core_set] = {}
    for year in [2013, 2014, 2015, 2016, 2017, 2018]:
        overlap_size_general[core_set][year] = {}
        overlap_size_refined[core_set][year] = {}
        overlap_general[core_set][year] = {}
        overlap_refined[core_set][year] = {}
        for cutoff in blast_clusters:
            overlap_general[core_set][year][cutoff] = core_blast_similar_pdbs[core_set][cutoff].intersection(data_sets[f'{year} general'])
            overlap_refined[core_set][year][cutoff] = core_blast_similar_pdbs[core_set][cutoff].intersection(data_sets[f'{year} refined'])
            overlap_size_general[core_set][year][cutoff] = len(overlap_general[core_set][year][cutoff]) - len(set(core_sets[core_set]).intersection(data_sets[f'{year} general']))
            overlap_size_refined[core_set][year][cutoff] = len(overlap_refined[core_set][year][cutoff]) - len(set(core_sets[core_set]).intersection(data_sets[f'{year} refined']))
    overlap_size_general[core_set] = pd.DataFrame.from_dict(overlap_size_general[core_set])
    overlap_size_refined[core_set] = pd.DataFrame.from_dict(overlap_size_refined[core_set])

data_set_overlap = {core_set: {data_set: {} for data_set in data_sets} for core_set in core_sets}
for core_set in core_sets:
    for cutoff in blast_clusters:
        for year in [2013, 2014, 2015, 2016, 2017, 2018]:
            data_set_overlap[core_set][f'{year} refined'][cutoff] = overlap_refined[core_set][year][cutoff]
            data_set_overlap[core_set][f'{year} general'][cutoff] = overlap_general[core_set][year][cutoff]

oob_rp_timesplit_seqid = {}
test_rp_timesplit_seqid = {}

# load list of structures whose ligand has a morgan fingerprint tanimoto coefficient >0.9 to the ligand of any structure in a core set
similar_ligands = {}
for year in ['2007', '2013', '2016']:
    with open(f'../data/pdbbind_2018_general_{year}_core_similar_ligands.txt') as f:
        similar_ligands[year] = [l.strip() for l in f]

similar_ligands['all'] = list(set(similar_ligands['2007'] + similar_ligands['2013'] + similar_ligands['2016']))

# Train with each feature set, using general and refined from each year, with and without structures containing similar ligands to those in the test set removed from the training set.
# Test on each core set, as well as the combined set of core structures
oob_rp = {}
test_rp = {}
oob_rp_ligands_removed = {}
test_rp_ligands_removed = {}
predictions = {}
predictions_ligands_removed = {}
# We also want to train with nothing removed from the training data
blast_clusters['None'] = []
for core_set in core_sets:
    for data_set in data_sets:
        data_set_overlap[core_set][data_set]['None'] = set()

n_fits = len(blast_clusters) * len(data_sets) * len(feature_sets) * len(core_sets)
i = 0

true_values = {}
test_sets = {}
for core_set in core_sets:
    test = [pdb for pdb in core_sets[core_set] if pdb in all_features.index]
    test_sets[core_set] = test
    true_values[core_set] = binding_data.loc[test].values.ravel().tolist()
with open('../results/test_set_binding_data.json', 'w') as f:
    json.dump(true_values, f)

with open('../data/core_set_pdbs.json', 'w') as f:
    json.dump(test_sets, f)

print('Fitting models')

for cutoff in blast_clusters:
    predictions[cutoff] = {}
    predictions_ligands_removed[cutoff] = {}
    oob_rp[cutoff] = {}
    test_rp[cutoff] = {}
    oob_rp_ligands_removed[cutoff] = {}
    test_rp_ligands_removed[cutoff] = {}

    for data_set in data_sets:
        predictions[cutoff][data_set] = {}
        predictions_ligands_removed[cutoff][data_set] = {}
        oob_rp[cutoff][data_set] = {}
        test_rp[cutoff][data_set] = {}
        oob_rp_ligands_removed[cutoff][data_set] = {}
        test_rp_ligands_removed[cutoff][data_set] = {}

        for feature_set in feature_sets:
            predictions[cutoff][data_set][feature_set] = {}
            predictions_ligands_removed[cutoff][data_set][feature_set] = {}
            oob_rp[cutoff][data_set][feature_set] = {}
            test_rp[cutoff][data_set][feature_set] = {}
            oob_rp_ligands_removed[cutoff][data_set][feature_set] = {}
            test_rp_ligands_removed[cutoff][data_set][feature_set] = {}

            for core_set in core_sets:

                # first without removing similar ligands
                overlap = data_set_overlap[core_set][data_set][cutoff]
                train = [pdb for pdb in all_features.index if pdb in data_sets[data_set] and pdb not in core_sets['all'] and pdb not in overlap]
                X_train = feature_sets[feature_set].loc[train]
                y_train = binding_data.loc[train]

                test = [pdb for pdb in core_sets[core_set] if pdb in all_features.index]
                X_test = feature_sets[feature_set].loc[test]
                y_test = binding_data.loc[test]

                model = RandomForestRegressor(n_estimators=500, max_features=0.33, n_jobs=64, random_state=42, oob_score=True)
                model.fit(X_train, y_train)
                filename = 'rf_all-ligs.sav'
                joblib.dump(model, filename)

                pred = list(model.predict(X_test))
                predictions[cutoff][data_set][feature_set][core_set] = pred
                
                oob_rp[cutoff][data_set][feature_set][core_set] = pearsonr(y_train.values.ravel(), model.oob_prediction_)[0]
                test_rp[cutoff][data_set][feature_set][core_set] = pearsonr(y_test.values.ravel(), pred)[0]

                # now remove similar ligands
                train = [pdb for pdb in train if pdb not in similar_ligands[core_set]]
                X_train = feature_sets[feature_set].loc[train]
                y_train = binding_data.loc[train]

                model = RandomForestRegressor(n_estimators=500, max_features=0.33, n_jobs=64, random_state=42, oob_score=True)
                model.fit(X_train, y_train)
                filename = 'rf_diverse-ligs.sav'
                joblib.dump(model, filename)

                pred = list(model.predict(X_test))
                predictions_ligands_removed[cutoff][data_set][feature_set][core_set] = pred

                oob_rp_ligands_removed[cutoff][data_set][feature_set][core_set] = pearsonr(y_train.values.ravel(), model.oob_prediction_)[0]
                test_rp_ligands_removed[cutoff][data_set][feature_set][core_set] = pearsonr(y_test.values.ravel(), pred)[0]

                i += 2
                if i % 50 == 0:
                    print(f'Fitted {i} of {2 * n_fits}.')

with open('../results/varying_seqid_cutoff_test_results.json', 'w') as f:
    json.dump(test_rp, f)
with open('../results/varying_seqid_cutoff_oob_results.json', 'w') as f:
    json.dump(oob_rp, f)
with open('../results/varying_seqid_cutoff_similar_ligands_removed_test_results.json', 'w') as f:
    json.dump(test_rp_ligands_removed, f)
with open('../results/varying_seqid_cutoff_similar_ligands_removed_oob_results.json', 'w') as f:
    json.dump(oob_rp_ligands_removed, f)

with open('../results/varying_seqid_cutoff_predictions.json', 'w') as f:
    json.dump(predictions, f)
with open('../results/varying_seqid_cutoff_similar_ligands_removed_predictions.json', 'w') as f:
    json.dump(predictions_ligands_removed, f)

