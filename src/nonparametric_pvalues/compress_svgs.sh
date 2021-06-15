#!/bin/bash
#SBATCH --job-name=compress_svg
#SBATCH --partition=dpwall
#SBATCH --output=/scratch/users/briannac/logs/compress_svg.out
#SBATCH --error=/scratch/users/briannac/logs/compress_svg.err
#SBATCH --time=40:00:00
#SBATCH --mem=75G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

### file at $MY_HOME/sequence_based_biomarkers/src/nonparametric_pvalues/compress_svg.sh

cd $MY_HOME
scour -i sequence_based_biomarkers/results/differential_abundances/obesity_clusters.svg -o sequence_based_biomarkers/results/differential_abundances/obesity_clusters_compressed.svg --enable-viewboxing --enable-id-stripping   --enable-comment-stripping --shorten-ids --indent=none


scour -i sequence_based_biomarkers/results/differential_abundances/autism_clusters.svg -o sequence_based_biomarkers/results/differential_abundances/autism_clusters_compressed.svg --enable-viewboxing --enable-id-stripping   --enable-comment-stripping --shorten-ids --indent=none
