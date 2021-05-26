#!/bin/sh
#SBATCH --job-name=compute_false_discoveries
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/users/briannac/logs/compute_false_discoveries_%a.out
#SBATCH --error=/scratch/users/briannac/logs/compute_false_discoveries_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=2TB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at /home/groups/dpwall/briannac/sequence_based_biomarkers/src/permutation_test/compute_false_discoveries.sh

module load python/3.6.1

for biomarker in asv taxa otu90 otu95 otu97 otu99 sbb1 sbb2 sbb3; do
    for dataset in autism obesity; do
        python3.6 -u $MY_HOME/sequence_based_biomarkers/src/permutation_test/compute_false_discoveries.py $dataset $biomarker
    done
done