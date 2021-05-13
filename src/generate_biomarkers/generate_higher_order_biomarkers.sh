#!/bin/bash
#SBATCH --job-name=generate_higher_order_biomarkers
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/users/briannac/logs/generate_higher_order_biomarkers.out
#SBATCH --error=/scratch/users/briannac/logs/generate_higher_order_biomarkers.err
#SBATCH --time=20:00:00
#SBATCH --mem=1000GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu


### file at $MY_HOME/sequence_based_biomarkers/src/generate_SBBs/generate_higher_order_biomarkers.sh
module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36
ml python/3.6
cd $MY_HOME/sequence_based_biomarkers

for study in yogurt; do #obese_lean_twins diabetes_150; do
    for msa_method in DECIPHER; do
        echo "removing duplicates"
        python3.6 src/generate_SBBs/remove_duplicate_variants.py $study ${msa_method}
        for order in 1; do #2 3; do
            echo running $study $msa_method $order
            echo "geneating higher order biomarkers"
            python3.6 src/generate_SBBs/sbbs.py $order data/${study}/person_vs_taxa.tsv data/${study}/${msa_method}_taxa_vs_variants_unique.tsv results/generate_SBBs/${study}_${msa_method}_
            echo "removing SBBs in LD"
            python3.6 -u src/generate_SBBs/remove_SBB_in_LD.py $study ${msa_method} $order
            \rm results/generate_SBBs/${study}_${msa_method}_person_variant${order}_condensed.npy
            \rm results/generate_SBBs/${study}_${msa_method}_biomarkers${order}.txt
            \rm results/generate_SBBs/${study}_${msa_method}_biomarker_exists${order}.npy
            \rm results/generate_SBBs/${study}_${msa_method}_biomarkers${order}_unique.npy
        done
    done
done


study=yogurt
msa_method=DECIPHER
order=1


