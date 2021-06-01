#!/bin/bash
#SBATCH --job-name=lefse
#SBATCH --partition=dpwall
#SBATCH --output=/scratch/users/briannac/logs/lefse.out
#SBATCH --error=/scratch/users/briannac/logs/lefse.err
#SBATCH --time=1:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/lefse/lefse.sh


ml R/3.6.1
ml python/3.6
for biomarker in sbb2 sbb3; do
    for dataset in autism obesity; do
    
        \rm $MY_HOME/sequence_based_biomarkers/results/lefse/sample_vs_biomarker_lefse_${biomarker}_${dataset}.out
        
        echo $biomarker $dataset
        python3.6 /oak/stanford/groups/dpwall/computeEnvironments/lefse/format_input.py \
        $MY_HOME/sequence_based_biomarkers/intermediate_files/lefse/sample_vs_biomarker_lefse_${biomarker}_${dataset}.tsv  $MY_HOME/sequence_based_biomarkers/intermediate_files/lefse/sample_vs_biomarker_lefse_${biomarker}_${dataset}.in \
        -c 1 -u 2  -f c 

        python3.6 /oak/stanford/groups/dpwall/computeEnvironments/lefse/run_lefse.py \
        $MY_HOME/sequence_based_biomarkers/intermediate_files/lefse/sample_vs_biomarker_lefse_${biomarker}_${dataset}.in \
        $MY_HOME/sequence_based_biomarkers/results/lefse/sample_vs_biomarker_lefse_${biomarker}_${dataset}.out -s 2 

    done
done
