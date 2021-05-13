#!/bin/sh
#SBATCH --job-name=compute_data_pvalues
#SBATCH --partition=owners
#SBATCH --array=1-7
#SBATCH --output=/scratch/users/briannac/logs/compute_data_pvalues_%a.out
#SBATCH --error=/scratch/users/briannac/logs/compute_data_pvalues_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at /home/groups/dpwall/briannac/sequence_based_biomarkers/src/permutation_test/compute_data_pvalues.sh

## SLURM_ARRAY_TASK_ID=1 ### CHANGE WHEN RUNNING BATCH!!!


module load python/3.6.1

N=$((SLURM_ARRAY_TASK_ID+i))
head -n $((N+1)) /home/groups/dpwall/briannac/sequence_based_biomarkers/intermediate_files/permuation_test/list_of_iters_unfinished.tsv | tail -n 1  > $SCRATCH/tmp/tmp_compute_data_pvalues_$SLURM_ARRAY_TASK_ID.txt

while read OUTFILE BIOMARKERS DATATYPE ; do
    python3.6 -u $MY_HOME/sequence_based_biomarkers/src/permutation_test/compute_data_pvalues.py \
    $BIOMARKERS $DATATYPE $OUTFILE
done < $SCRATCH/tmp/tmp_compute_data_pvalues_$SLURM_ARRAY_TASK_ID.txt
\rm $SCRATCH/tmp/tmp_compute_data_pvalues_$SLURM_ARRAY_TASK_ID.txt
