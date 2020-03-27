# Learning ligand
This repository was hacked together from the work in "Learning from the ligand: using ligand-based features to improve binding affinity prediction" (https://academic.oup.com/bioinformatics/article-abstract/36/3/758/5554651). The starting point was download of the data and code found at http://opig.stats.ox.ac.uk/resources (as referenced in the manuscript)

The added scripts were used specifically to a) attempt to re-build the models and b) to use the rebuilt models to score a number of different ligand poses against various protein structures.

In the very near future I will clean this up so that the whole pipeline is runnable in an easy number of short steps 

NB: the generated model files and results for our use-case are not here. They are too big. 

## Re-building (roughly) the models from the paper
The script in `scripts/rebuild_models.py` should do this. It takes the feature sets that were distributed on the space mentioned in the paper and re-runs sci-kit learn to build the various models. These are output in the `models` directory as `.sav` files

## Running the models on your own data

### Cleaning data and calculating features
For our use-case, we had many ligands to be scored against one protein, and this case many times. The script `scripts/process-mpro.py` shows an example of how to do this. The script:
1. Cleans up the input file (very specific to our case)
2. Splits the input file by the protein target all of the ligands (again, quite specific)
3. Calculates the rdkit, rfscore and binana features that are required for scoring with the models we have re-built

### Predicting PkA for ligands with models we re-built
This is again quite sppecific for our use case, but an example is shown in `scripts/run_models_mpro.py`. The script:
1. loads the feature sets
2. calculates the scores for each set of ligands
3. writes the scores to the input sdf file
