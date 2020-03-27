"""Computes BINANA features for the PDBbind data set.

Usage:
    compute_binana_features.py [-vh] <pdb_list_file> <pdbbind_dir> <output_file> [<blacklist_file>] [--num-cores=<n>]

Arguments:
    pdb_list_file   file containing pdb codes of complexes to use
    pdbbind_dir     top-level directory of the PDBbind 2017 database
    output_file     file to save the computed features to
    blacklist_file  optional file containing PDB codes to skip; some PDB
        structures are known to cause OpenBabel to seg fault

Options:
    -h --help             show this message and exit
    -v --verbose          print progress updates
    --num-cores=<n>  number of cores to use [default: -1]

Computes all BINANA features as implemented in Open Drug Discovery Toolkit
for specified complexes from the PDBbind data set and saves the results to
the specified file in .csv format.

"""
import pathlib
import glob
import pandas as pd

from joblib import delayed, Parallel
from oddt.toolkits import ob
from oddt.scoring.descriptors import binana


def featurise(protein_file, ligand_file):
    """Compute BINANA features for a protein and a ligand and return the results as a dictionary.

    Args:
        protein_file (str): Name of the .pdb file containing the protein.
        ligand_file (str): Name of the .sdf file containing the ligand.

    Returns:
        result (dict): A dictionary containing the BINANA features for the protein-ligand complex.
    """

    protein = next(ob.readfile('pdb', protein_file))
    protein.protein = True # DO NOT SKIP THIS or you will be a sad panda
    ligand = next(ob.readfile('sdf', ligand_file))

    # Build BINANA features
    binana_engine = binana.binana_descriptor(protein)
    result = {name: value for name, value in zip(binana_engine.titles, binana_engine.build([ligand])[0])}
    return result

def binana_multi_ligs(protein_file, ligands_file):
    protein = next(ob.readfile('pdb', protein_file))
    protein.protein = True  # DO NOT SKIP THIS or you will be a sad panda
    ligs = list(ob.readfile('sdf', ligands_file))

    fname = ligands_file.split('/')[-1].replace('.sdf', '')

    results = {}

    for i, ligand in enumerate(ligs):
        binana_engine = binana.binana_descriptor(protein)
        result = {name: value for name, value in zip(binana_engine.titles, binana_engine.build([ligand])[0])}
        results[str(fname + '_' + str(i))] = result

    return results



if __name__=='__main__':
    args = docopt(__doc__)

    pdb_list_file = args['<pdb_list_file>']
    pdbbind_dir = args['<pdbbind_dir>']
    output_file = args['<output_file>']
    blacklist_file = args['<blacklist_file>']

    verbose = 10 if args['--verbose'] else 0
    n_jobs = int(args['--num-cores'])

    with open(pdb_list_file, 'r') as f:
        pdbs = [l.strip() for l in f]

    # drop any blacklisted PDB codes
    if blacklist_file:
        with open(os.path.join(blacklist_file), 'r') as f:
            blacklist = [line.strip() for line in f]
        pdbs = [pdb for pdb in pdbs if pdb not in blacklist]

    # list protein and ligand pdb/sdf files and compute features
    protein_files = {pdb: str(pathlib.Path(pdbbind_dir, pdb, f'{pdb}_protein.pdb')) for pdb in pdbs}
    ligand_files = {pdb: str(pathlib.Path(pdbbind_dir, pdb, f'{pdb}_ligand.sdf')) for pdb in pdbs}
    with Parallel(n_jobs=n_jobs, verbose=verbose) as parallel:
        results = parallel(delayed(featurise)(protein_files[pdb], ligand_files[pdb]) for pdb in pdbs)

    features = {pdb: result for pdb, result in zip(pdbs, results)}
    features = pd.DataFrame.from_dict(features).T
    features.to_csv(output_file)

