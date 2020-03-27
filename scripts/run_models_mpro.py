import pandas as pd
import glob
import joblib
import os
from rdkit import Chem

models_rdkit_features_cols = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_rdkit_features_clean.csv'), index_col=0).columns
models_rfscore_features_cols = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_rfscore_features_clean.csv'), index_col=0).columns
models_nnscore_features_cols = pd.read_csv(os.path.join('..', 'data', 'pdbbind_2018_general_binana_features_clean.csv'), index_col=0).columns

for f in glob.glob('../files_to_analyse/Mpro-x*_mols.sdf'):
    if os.path.isfile(f.replace('_mols', '_mols_scored')):
        continue
    fname = f.split('/')[-1].replace('.sdf', '').replace('_mols', '')
    rfscore_features = pd.read_csv(f'../files_to_analyse/features/rfscore-features_{fname}.csv', index_col=0)
    rfscore_features = rfscore_features[models_rfscore_features_cols]
    rdkit_features = pd.read_csv(f'../files_to_analyse/features/rdkit-features_{fname}.csv', index_col=0)
    rdkit_features = rdkit_features[models_rdkit_features_cols]
    nnscore_features = pd.read_csv(f'../files_to_analyse/features/binana-features_{fname}.csv', index_col=0)
    nnscore_features = nnscore_features[models_nnscore_features_cols]

    all_features = pd.concat([rdkit_features, rfscore_features, nnscore_features], axis='columns')

    print(nnscore_features.columns)
    print(models_nnscore_features_cols)

    names = list(all_features.index)

    feature_sets = {
        'RDKit': rdkit_features,
        'RF-Score': rfscore_features,
        'RF-Score + RDKit': pd.concat([rdkit_features, rfscore_features], axis='columns'),
        'NNScore 2.0': nnscore_features,
        'NNScore 2.0 + RDKit': pd.concat([rdkit_features, nnscore_features], axis='columns')

    }

    vina_terms = pd.Index(
        ['vina_gauss1', 'vina_gauss2', 'vina_hydrogen', 'vina_hydrophobic', 'vina_repulsion', 'num_rotors'])
    feature_sets['Vina'] = all_features[vina_terms]
    feature_sets['Vina + RDKit'] = all_features[vina_terms.union(rdkit_features.columns)]
    feature_sets['Vina + RDKit'] = feature_sets['Vina + RDKit'].drop('NumRotatableBonds', axis='columns')
    feature_sets['RF-Score v3'] = all_features[vina_terms.union(rfscore_features.columns)]
    feature_sets['RF-Score v3 + RDKit'] = pd.concat([feature_sets['RF-Score v3'], feature_sets['RDKit']],
                                                    axis='columns')

    print(feature_sets['RF-Score v3 + RDKit'])

    scores = {}

    # for model in models:
    for feature_set in feature_sets:
        models = glob.glob(f'../models/*set_{feature_set}.sav')
        for model in models:
            print(model)
            score_name = ''.join(model.split('/')[-1].rsplit()).replace('+', '-').replace('.sav', '').replace('model', 'score')
            print(score_name)
            print(feature_set)
            loaded_model = joblib.load(model)
            try:
                X_test = [feature_sets[feature_set].loc[name] for name in names]
                result = loaded_model.predict(X_test)
            except ValueError:
                print(f'failed {model} {fname}')
                break

            score_iter = list(zip(names, result))

            # print(score_iter)
            # print(result)
            scores[score_name] = score_iter

    print(len(scores))

    sdf_file = f'../files_to_analyse/{fname}_mols.sdf'
    suppl = Chem.SDMolSupplier(sdf_file)
    mols = [m for m in suppl]

    for k in scores.keys():
        for mol, val in scores[k]:
            i = mol.split('_')[-1]
            print(val)
            mols[int(i)].SetProp(k, str(round(val, 3)))

    out_sdf = f'../files_to_analyse/{fname}_mols_scored.sdf'
    w = Chem.SDWriter(out_sdf)
    for mol in mols:
        w.write(mol)





