# Generate SBBs


## Storage
**Intermediate files**: ``` ```
**Archived Intermediate files**: ``` ```
**Final**: ``` ```

## Prerequisites
0.1. Run preprocessing pipeline.

0.2 

## Pipeline 

Outputs to be 
- ```sample_vs_biomarker_<BIOMARKER>_<DISEASE>.tsv```
- ```taxa_vs_biomarker_<BIOMARKER>_<DISEASE>.tsv```


1. ```align.sh```
    - ***Inputs***: ```seqs.fa```
    - ***Outputs***: ```aligned_seqs.fa```

2. ```format_data.ipynb```:
    - ***Inputs***: ```aligned_seqs.fa```
    - ***Outputs***: ```otu_snp_table.tsv```

3. ```sbbs.sh```:
    - ***Inputs***: ```otu_snp_table.tsv```, ```sample_vs_otu.tsv```
    - ***Outputs***: ```sample_vs_sbb<order>.tsv```, ```otu_vs_sbb<order>.tsv
    
4. ```comparison_methods.ipynb```: 
    - ***Inputs***: ```otu_snp_table.tsv```, ```sample_vs_otu.tsv```
    - ***Outputs***: ```sample_vs_sbb<order>_unique.tsv```, ```otu_vs_sbb<order>_unique.tsv
