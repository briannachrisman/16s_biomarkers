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

for biomarker in sbb3; do
    for dataset in autism; do
        echo $dataset $biomarker
        python3.6 $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/rank_stat.py $dataset $biomarker
    done
done