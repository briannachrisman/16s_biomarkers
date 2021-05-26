#!/bin/bash
#SBATCH --job-name=sbbs_autism
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/users/briannac/logs/sbbs_autism3.out
#SBATCH --error=/scratch/users/briannac/logs/sbbs_autism3.err
#SBATCH --time=20:00:00
#SBATCH --mem=2T
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu


### file at $MY_HOME/sequence_based_biomarkers/src/generate_biomarkers/sbbs.sh
module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36
ml python/3.6

for order in 2 3; do
    for dataset in autism_small obesity_small; do
        echo $dataset $order SBBs
        python3.6 -u $MY_HOME/sequence_based_biomarkers/src/generate_biomarkers/sbbs.py $order $dataset
    done
done