#!/bin/bash
#SBATCH --job-name=lefse_format
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/users/briannac/logs/lefse_format.out
#SBATCH --error=/scratch/users/briannac/logs/lefse_format.err
#SBATCH --time=40:00:00
#SBATCH --mem=2TB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/phenotype_associations/lefse_format.sh

python3.6 $MY_HOME/sequence_based_biomarkers/src/phenotype_associations/lefse_format.py