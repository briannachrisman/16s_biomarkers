#!/bin/bash
#SBATCH --job-name=fdr
#SBATCH --partition=dpwall
#SBATCH --array=1-10%5
#SBATCH --output=/scratch/users/briannac/logs/fdr.out
#SBATCH --error=/scratch/users/briannac/logs/fdr.err
#SBATCH --time=40:00:00
#SBATCH --mem=75G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/fdr.sh

ml python/3.6
#START=$((SLURM_ARRAY_TASK_ID*100-100))
#END=$((START+100))
for biomarker in otu95 otu97 otu99 micropheno4 micropheno6 micropheno8 sbb1 sbb2 sbb3; do #  ; do
    for dataset in obesity autism; do
        echo $dataset $biomarker
        mkdir $MY_HOME/sequence_based_biomarkers/intermediate_files/nonparametric_pvalues/${dataset}_${biomarker}
        echo "rank stat..."
        python3.6 $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/rank_stat.py $dataset $biomarker
        echo "permutation test..."
        python3.6 -u $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/fdr.py $dataset $biomarker $SLURM_ARRAY_TASK_ID
    done
done
