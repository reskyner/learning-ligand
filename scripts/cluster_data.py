"""Clusters the PDBbind 2017 refined set by 90% similarity to proteins in the core set.

Usage:
    cluster_data.py [-h] <pdbbind_2016_dir> <cluster_file> <training_set_file> <test_set_file>

Arguments:
    pdbbind_2016_dir     top-level directory of the PDBbind 2016 database
    cluster_file         file containing BLASTClust clusters of the PDBbind 2017 refined set in .json format
    training_set_file    file containing PDB codes of the training set
    test_set_file        file containing PDB codes of the test set

Options:
    -h --help       show this message and exit

This script is gnarly and inefficient, but documented and functional.
It is neither well-designed nor well-engineered. I'm sorry.
Embrace the chaos of throwaway research code. Climb the ladder.

It is assumed that the PDBbind 2016 indices live in the standard directory structure
i.e. <top-level-dir>/index/<filename>. If you just downloaded the core cluster file
and saved it somewhere else, go ahead and edit the script. You'll probably make it better.

Clusters of the refined set corresponding to the 56 clusters of the core set we use
are saved to ../data/single_target_clusters.json

A list of PDB codes belonging to the training set whose proteins are in the same
BLASTClust clusters as any test set protein is written to ../data/training_test_overlap.txt

"""
import json
import os

from docopt import docopt


# parse command line arguments
args = docopt(__doc__)
pdbbind_2016_dir = args['<pdbbind_2016_dir>']
cluster_file = args['<cluster_file>']
test_set_file = args['<test_set_file>']
training_set_file = args['<training_set_file>']

# parse protein names and clustering data from PDBbind indices

core_cluster_file = os.path.join(pdbbind_2016_dir, 'index', 'INDEX_core_cluster.2016')
with open(core_cluster_file, 'r') as f:
    lines = [l.strip().split()[0] for l in f if not l.startswith('#')]
core_cluster_data = [lines[x:x+5] for x in range(0, len(lines),5)]

core_names_file = os.path.join(pdbbind_2016_dir, 'index', 'INDEX_core_name.2016')
with open(core_names_file, 'r') as f:
    lines = [line.strip().split() for line in f if not line.startswith('#')]
core_protein_names = {line[0]: ' '.join(line[3:]) for line in lines}

# dict of protein name: pdb codes for each cluster of the core set
# note - there are actually two clusters for beta-lactamase
core_cluster_names = list(set([name for _, name in core_protein_names.items()]))
core_clusters = {name: [] for name in core_cluster_names}
for pdb in core_protein_names:
    core_clusters[core_protein_names[pdb]].append(pdb)

# training and test set PDB codes
with open(training_set_file, 'r') as f:
    training_set = [line.strip() for line in f]

with open(test_set_file, 'r') as f:
    test_set = [line.strip() for line in f]

# clusters of the full data set
with open(cluster_file, 'r') as f:
    clusters = json.load(f)

# using sets speeds up checking membership of each cluster
set_clusters = [set(c) for c in clusters]

# for each core set cluster, identify all proteins in the refined set which
# are >90% sequence identical to any protein from the core set cluster
refined_clusters = {}
for c in core_clusters:

    # first, add the proteins from the core set that are in our test set
    refined_clusters[c] = [i for i in core_clusters[c] if i in test_set]

    # for each protein from the core set, add all proteins from each refined
    # set cluster to which the core set protein belongs
    for pdb in core_clusters[c]:
        for cluster in set_clusters:
            if pdb in cluster:
                refined_clusters[c].extend([i for i in cluster if i in training_set])

    # remove repeated entries
    refined_clusters[c] = list(set(refined_clusters[c]))

# finally, list the PDB codes of all proteins in the training data that
# belong to the same cluster as any protein in the test set
similar_pdbs = [pdb for c in refined_clusters for pdb in refined_clusters[c]]
training_set_overlap = list(set(similar_pdbs).difference(set(test_set)))

# finally, dump all the clustering information and never look back

with open(os.path.join('..', 'data', 'training_set_overlap.txt'), 'w') as f:
    print(*training_set_overlap, sep='\n', file=f)

with open(os.path.join('..', 'data', 'single_target_clusters.json'), 'w') as f:
    json.dump(refined_clusters, f)

