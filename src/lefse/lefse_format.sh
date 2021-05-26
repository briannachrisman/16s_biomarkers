#!/bin/bash
#SBATCH --job-name=lefse_format
#SBATCH --partition=dpwall
#SBATCH --output=/scratch/users/briannac/logs/lefse_format.out
#SBATCH --error=/scratch/users/briannac/logs/lefse_format.err
#SBATCH --time=20:00:00
#SBATCH --mem=128GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/lefse/lefse_format.sh

python3.6 -u $MY_HOME/sequence_based_biomarkers/src/lefse/lefse_format.py