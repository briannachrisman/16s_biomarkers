#!/bin/bash
#SBATCH --job-name=rank_stat_autism
#SBATCH --partition=dpwall
#SBATCH --output=/scratch/users/briannac/logs/rank_stat_autism.out
#SBATCH --error=/scratch/users/briannac/logs/rank_stat_autism.err
#SBATCH --time=40:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/rank_stat_autism.sh

for biomarker in taxa asv otu90 otu95 otu97 otu99 micropheno4 micropheno6 micropheno8 sbb1 sbb2 sbb3; do
    for dataset in autism; do
        echo $dataset $biomarker
        python3.6 -u $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/rank_stat.py $dataset $biomarker
    done
done