#!/bin/bash
#SBATCH --job-name=RFC
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/users/briannac/logs/RFC.out
#SBATCH --error=/scratch/users/briannac/logs/RFC.err
#SBATCH --time=3:00:00
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/ML_phenotype_prediction/logit.sh

for biomarker in taxa asv otu90 otu95 otu97 otu99 micropheno4 micropheno6 micropheno8 sbb1 sbb2 sbb3; do
    for dataset in autism obesity; do
        echo $dataset $biomarker
        python3.6 -u $MY_HOME/sequence_based_biomarkers/src/ML_phenotype_prediction/RFC.py $dataset $biomarker
    done
done
