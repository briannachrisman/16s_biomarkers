#!/bin/bash
#SBATCH --job-name=align
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/users/briannac/logs/align.out
#SBATCH --error=/scratch/users/briannac/logs/align.err
#SBATCH --time=20:00:00
#SBATCH --mem=1000GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu


### file at $MY_HOME/sequence_based_biomarkers/src/generate_biomarkers/align.sh

echo "Aligning...."
augur align \
  --sequences /home/groups/dpwall/briannac/sequence_based_biomarkers/data/autism/seqs.fa \
  --output /home/groups/dpwall/briannac/sequence_based_biomarkers/results/generate_biomarkers/seqs_aligned_autism_mafft.fa \
  --nthreads=auto
  
  
  
  echo "Aligning...."
augur align \
  --sequences /home/groups/dpwall/briannac/sequence_based_biomarkers/data/obesity/seqs.fa \
  --output /home/groups/dpwall/briannac/sequence_based_biomarkers/results/generate_biomarkers/seqs_aligned_obesity_mafft.fa \
  --nthreads=auto