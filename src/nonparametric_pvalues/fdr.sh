#!/bin/bash
#SBATCH --job-name=fdr
#SBATCH --partition=dpwall
#SBATCH --array=1
#SBATCH --output=/scratch/users/briannac/logs/fdr.out
#SBATCH --error=/scratch/users/briannac/logs/fdr.err
#SBATCH --time=40:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at 


START=$((SLURM_ARRAY_TASK_ID*100-100))
END=$((START+100))
for biomarker in sbb2 sbb3; do #  ; do
    for dataset in autism obesity; do
        echo $dataset $biomarker
        mkdir $MY_HOME/sequence_based_biomarkers/intermediate_files/nonparametric_pvalues/${dataset}_${biomarker}
        python3.6 -u $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/fdr.py $dataset $biomarker $START $END
    done
done