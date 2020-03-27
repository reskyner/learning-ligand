"""Computes RDKit features for the PDBbind data set.

Usage:
    compute_rdkit_features.py [-h] <pdb_list_file> <pdbbind_dir> <output_file>

Arguments:
    pdb_list_file   file containing pdb codes of complexes to use
    pdbbind_dir     top-level directory of the PDBbind data set
    output_file     file to save the computed features to

Options:
    -h --help       show this message and exit

Computes all 2D molecular descriptors as implemented in the Descrptors module
of RDKit for the ligands of the specified complexes from the PDBbind data set
and saves the results to the specified file in .csv format.

"""
import pathlib
import pandas as pd
import glob

# from docopt import docopt
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Descriptors

def fergus_rdkit():
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    args = docopt(__doc__)
    pdb_list_file = args['<pdb_list_file>']
    pdbbind_dir = args['<pdbbind_dir>']
    output_file = args['<output_file>']

    with open(pdb_list_file, 'r') as f:
        pdbs = [l.strip() for l in f]

    # load ligands and compute features
    features = {}
    descriptors = {d[0]: d[1] for d in Descriptors.descList}
    for pdb in pdbs:
        # prefer to use the .sdf provided by PDBbind
        sdf = pathlib.Path(pdbbind_dir, pdb, f'{pdb}_ligand.sdf')
        mol = next(Chem.SDMolSupplier(str(sdf), removeHs=False))

        # but we'll try the .mol2 if RDKit can't parse the .sdf
        if mol is None:
            mol2 = pathlib.Path(pdbbind_dir, pdb, f'{pdb}_ligand.mol2')
            mol = Chem.MolFromMol2File(str(mol2), removeHs=False)

        # skip the ligand if RDKit can't parse the .mol2
        if mol is None:
            continue

        try:
            features[pdb] = {d: descriptors[d](mol) for d in descriptors}
        except ValueError as e:
            print(e)
            continue

    features = pd.DataFrame.from_dict(features).T
    features.to_csv(output_file)

def rdkit_features_multilig(sdf):
    descriptors = {d[0]: d[1] for d in Descriptors.descList}
    features={}
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mols = [mol for mol in suppl]
    for i, mol in enumerate(mols):
        name = sdf.split('/')[-1].replace('.sdf', '') + '_' + str(i)
        features[name] = {d: descriptors[d](mol) for d in descriptors}
    return features



