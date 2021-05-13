# Sequence based biomarkers

## Prerequisites


## Pipeline 

1. ✓ ```generate_biomarkers.sh```:
    - **Inputs**: 
    - **Outputs**: ```results/generate_sbbs/<DATASET>_<MSA>_biomarkers<ORDER>_unique.npy```


2.   ```compute_data_p_values.sh```: 
    - **Inputs**:```results/generate_sbbs/<DATASET>_<MSA>_biomarkers<ORDER>_unique.npy```, ```data/<DATASET>/sample_data.tsv```
    - **Outputs**: ```results/permutation_test/<DATASET>_<MSA>_sbb<ORDER>_pvalues.tsv```


3.   ```permute_data.py```: Creates a permuted dataset.
    - **Inputs**:```results/generate_sbbs/<DATASET>_<MSA>_biomarkers<ORDER>_unique.npy```, ```data/<DATASET>/sample_data.tsv```
    - **Outputs**: ```intermediate_files/permutation_test/<DATASET>_permuted_phenos.tsv```


4.   ```compute_false_discoveries.sh```: 
    - **Inputs**: ```intermediate_files/covWAS/<CHROMOSOME>.<IMPROPER|PROPER|UNMAPPED>.tsv.gz>```
    - **Outputs**:  ```results/permutation_test/<DATASET>_<MSA>_sbb<ORDER>_pvalues.tsv```, ```intermediate_files/permutation_test/<DATASET>_permuted_phenos.tsv```, ```results/covWAS/<CHROMOSOME>.<IMPROPER|PROPER|UNMAPPED>.pvals_manhattan.svg>```
    
5.   ```concat_false_discoveries.sh```: 
    - **Inputs**: ```intermediate_files/covWAS/<CHROMOSOME>.<IMPROPER|PROPER|UNMAPPED>.tsv.gz>```
    - **Outputs**:  ```results/permutation_test/<DATASET>_<MSA>_sbb<ORDER>_pvalues.tsv```, ```intermediate_files/permutation_test/<DATASET>_permuted_phenos.tsv```, ```results/covWAS/<CHROMOSOME>.<IMPROPER|PROPER|UNMAPPED>.pvals_manhattan.svg>```
    
    
    
    