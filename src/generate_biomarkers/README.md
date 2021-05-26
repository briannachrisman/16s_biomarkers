# Generate SBBs


## Storage
**Intermediate files**: ``` ```
**Archived Intermediate files**: ``` ```
**Final**: ``` ```

## Prerequisites
0.1. Run preprocessing pipeline.

## Pipeline 

1. ```sbbs.sh```:
    - ***Inputs***: ```otu_snp_table.tsv```, ```sample_vs_otu.tsv```
    - ***Outputs***: ```sample_vs_sbb<order>.tsv```, ```otu_vs_sbb<order>.tsv
    
2. ```comparison_methods.ipynb```: 
    - ***Inputs***: ```otu_snp_table.tsv```, ```sample_vs_otu.tsv```
    - ***Outputs***: ```sample_vs_sbb<order>_unique.tsv```, ```otu_vs_sbb<order>_unique.tsv