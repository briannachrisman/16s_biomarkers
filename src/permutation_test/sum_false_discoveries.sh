#!/bin/sh
#SBATCH --job-name=compute_p_value_1sbb_autism
#SBATCH --partition=owners
#SBATCH --array=1-7
#SBATCH --output=/scratch/users/briannac/logs/compute_p_value_1sbb_autism_%a.out
#SBATCH --error=/scratch/users/briannac/logs/compute_p_value_1sbb_autism_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at /home/groups/dpwall/briannac/sequence_based_biomarkers/src/permutation_test/compute_p_value_1sbb_autism.sh

## SLURM_ARRAY_TASK_ID=1 ### CHANGE WHEN RUNNING BATCH!!!

cd $MY_HOME/sequence_based_biomarkers

cat intermediate_files/permutation_test/autism_sbb1_*.txt > intermediate_files/permutation_test/autism_sbb1.tsv
cat intermediate_files/permutation_test/autism_sbb2_*.txt > intermediate_files/permutation_test/autism_sbb2.tsv
cat intermediate_files/permutation_test/autism_sbb3_*.txt > intermediate_files/permutation_test/autism_sbb3.tsv
cat intermediate_files/permutation_test/autism_ASV_*.txt > intermediate_files/permutation_test/autism_ASV.tsv
cat intermediate_files/permutation_test/autism_ditaxa_*.txt > intermediate_files/permutation_test/autism_ditaxa.tsv
cat intermediate_files/permutation_test/autism_tax_levels_*.txt > intermediate_files/permutation_test/autism_tax_levels.tsv
cat intermediate_files/permutation_test/obesity_sbb1_*.txt > intermediate_files/permutation_test/obesity_sbb1.tsv
cat intermediate_files/permutation_test/obesity_sbb2_*.txt > intermediate_files/permutation_test/obesity_sbb2.tsv
cat intermediate_files/permutation_test/obesity_sbb3_*.txt > intermediate_files/permutation_test/obesity_sbb3.tsv
cat intermediate_files/permutation_test/obesity_ASV_*.txt > intermediate_files/permutation_test/obesity_ASV.tsv
cat intermediate_files/permutation_test/obesity_ditaxa_*.txt > intermediate_files/permutation_test/obesity_ditaxa.tsv
cat intermediate_files/permutation_test/obesity_tax_levels_*.txt > intermediate_files/permutation_test/obesity_tax_levels.tsv