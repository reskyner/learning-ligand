"""
Quick and dirty clustering of the PDBbind 2018 general/refined sets at 90% sequence identity
with respect to the 2007, 2013, and 2016 core sets
"""
import json
import pathlib

with open(pathlib.Path('..', 'data', 'pdb_blastclust90_clusters.json')) as f:
    pdb_clusters = json.load(f)

pdb_clusters = [set([pdb.lower() for pdb in c]) for c in pdb_clusters]

with open(pathlib.Path('..', 'data', 'pdbbind_2007_core_pdbs.txt')) as f:
    core_2007 = [l.strip() for l in f]

with open(pathlib.Path('..', 'data', 'pdbbind_2013_core_pdbs.txt')) as f:
    core_2013 = [l.strip() for l in f]

with open(pathlib.Path('..', 'data', 'pdbbind_2016_core_pdbs.txt')) as f:
    core_2016 = [l.strip() for l in f]

with open(pathlib.Path('..', 'data', 'pdbbind_2018_general_pdbs.txt')) as f:
    pdbbind = [l.strip() for l in f]

with open(pathlib.Path('..', 'data', 'pdbbind_2018_refined_pdbs.txt')) as f:
    refined = [l.strip() for l in f]
core_2007_similar_pdbs = []
core_2013_similar_pdbs = []
core_2016_similar_pdbs = []
for pdb in core_2007:
    for cluster in pdb_clusters:
        if pdb in cluster:
            core_2007_similar_pdbs.extend(cluster)
core_2007_similar_pdbs = sorted(list(set(core_2007_similar_pdbs)))

for pdb in core_2013:
    for cluster in pdb_clusters:
        if pdb in cluster:
            core_2013_similar_pdbs.extend(cluster)
core_2013_similar_pdbs = sorted(list(set(core_2013_similar_pdbs)))

for pdb in core_2016:
    for cluster in pdb_clusters:
        if pdb in cluster:
            core_2016_similar_pdbs.extend(cluster)
core_2016_similar_pdbs = sorted(list(set(core_2016_similar_pdbs)))

core_all_similar_pdbs = sorted(list(set(core_2007_similar_pdbs + core_2013_similar_pdbs + core_2016_similar_pdbs)))

core_all = sorted(list(set(core_2007 + core_2013 + core_2016)))

with open(pathlib.Path('..', 'data', 'pdbbind_2007_core_blastclust90_overlap.txt'), 'w') as f:
    print(*core_2007_similar_pdbs, sep='\n', file=f)

with open(pathlib.Path('..', 'data', 'pdbbind_2013_core_blastclust90_overlap.txt'), 'w') as f:
    print(*core_2013_similar_pdbs, sep='\n', file=f)

with open(pathlib.Path('..', 'data', 'pdbbind_2016_core_blastclust90_overlap.txt'), 'w') as f:
    print(*core_2016_similar_pdbs, sep='\n', file=f)

refined_overlap = [pdb for pdb in refined if pdb in core_2007_similar_pdbs and pdb not in core_2007]
print(len(refined_overlap))
refined_overlap = [pdb for pdb in refined if pdb in core_2013_similar_pdbs and pdb not in core_2013]
print(len(refined_overlap))
refined_overlap = [pdb for pdb in refined if pdb in core_2016_similar_pdbs and pdb not in core_2016]
print(len(refined_overlap))

general_overlap = [pdb for pdb in pdbbind if pdb in core_2007_similar_pdbs and pdb not in core_2007]
print(len(general_overlap))
general_overlap = [pdb for pdb in pdbbind if pdb in core_2013_similar_pdbs and pdb not in core_2013]
print(len(general_overlap))
general_overlap = [pdb for pdb in pdbbind if pdb in core_2016_similar_pdbs and pdb not in core_2016]
print(len(general_overlap))

refined_overlap = [pdb for pdb in refined if pdb in core_all_similar_pdbs and pdb not in core_all]
print(len(refined_overlap))
general_overlap = [pdb for pdb in pdbbind if pdb in core_all_similar_pdbs and pdb not in core_all]
print(len(general_overlap))

