"""Computes Tanimoto similarity between ligands in PDBbind general and core sets

"""
import json
import pathlib

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint

mols = {}
fingerprints = {}
test_pdbs = {}

for l in test_ligands:
    pdbs = [pdb for pdb, _ in ligand_groups[l] if pdb in binding_data.index]
    for pdb in pdbs:
        src = f'/data/griffin/fboyles/pdbbind_2017/{pdb}/{pdb}_ligand.sdf'
        mol = next(Chem.SDMolSupplier(src))
        if mol is None:
            src = f'/data/griffin/fboyles/pdbbind_2017/{pdb}/{pdb}_ligand.mol2'
            mol = Chem.MolFromMol2File(src)
        if mol:
            mols[l] = mol
            fingerprints[l] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            test_pdbs[l] = pdb
            break

all_mols = {}
all_fps = {}

for pdb in binding_data.index:
    src = f'/data/griffin/fboyles/pdbbind_2017/{pdb}/{pdb}_ligand.sdf'
    mol = next(Chem.SDMolSupplier(src))
    if mol is None:
        src = f'/data/griffin/fboyles/pdbbind_2017/{pdb}/{pdb}_ligand.mol2'
        mol = Chem.MolFromMol2File(src)
    all_mols[pdb] = mol
    all_fps[pdb] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

similarity = {}

for pdb1 in binding_data.index:
    similarity[pdb1] = {}
    for pdb2 in binding_data.index:
        similarity[pdb1][pdb2] = DataStructs.FingerprintSimilarity(all_fps[pdb1], all_fps[pdb2])
        
mapper = {test_pdbs[l]: l for l in test_pdbs}
similarity_df = pd.DataFrame.from_dict(similarity).loc[:,[test_pdbs[l] for l in test_pdbs]].rename(mapper=mapper, axis='columns')

similarity_df.to_csv('pdbbind_2017_general_ligand_tc.csv')
with open('lbap_test_pdbs.json', 'w') as f:
    json.dump(test_pdbs, f)
