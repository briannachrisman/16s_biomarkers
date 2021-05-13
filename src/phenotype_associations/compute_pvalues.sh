#!/bin/bash
#SBATCH --job-name=compute_pvalues
#SBATCH --partition=dpwall
#SBATCH --output=/scratch/users/briannac/logs/compute_pvalues.out
#SBATCH --error=/scratch/users/briannac/logs/compute_pvalues.err
#SBATCH --time=40:00:00
#SBATCH --mem=128GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/phenotype_associations/compute_pvalues.sh

for biomarker in sbb1 sbb2 sbb3 taxa asv otu90 otu95 otu97 otu99 micropheno; do
    for dataset in autism obesity; do
        python3.6 $MY_HOME/sequence_based_biomarkers/src/phenotype_associations/compute_pvalues.py $dataset $biomarker
    done
done