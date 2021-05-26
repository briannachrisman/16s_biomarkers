# File at $MY_HOME/sequence_based_biomarkers/src/permutation_test/generate_permuted_data.py

# Set up libraries, etc.
import pandas as pd
import numpy as np
import seaborn as sns
from collections import Counter
import scipy.stats as ss
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
from scipy import stats
PROJECT_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'

#obesity_person_vs_taxa = pd.read_table(PROJECT_DIR + 'data/obesity/sample_vs_asv.tsv').transpose()
obesity_metadata =  pd.read_table(PROJECT_DIR + 'data/obesity/sample_metadata.tsv', index_col=0)

#autism_person_vs_taxa = pd.read_table(PROJECT_DIR + 'data/autism/sample_vs_asv.tsv').transpose()
autism_metadata =  pd.read_table(PROJECT_DIR + 'data/autism/sample_metadata.tsv', index_col=0)

np.random.seed(42)

def Generate_Random_Table(metadata, n_iters=1000, datatype='autism'):
    phenos_table = np.zeros((n_iters, len(metadata)))
    for n_iter in tqdm.tqdm(range(n_iters)):
        if datatype=='autism':
            swap_dict = {f:np.random.random()>.5 for f in np.unique(metadata.family)}
            new_phenos = [abs(phenotype-swap_dict[family]) for family, phenotype in zip(metadata.family, metadata.phenotype)]
            while sum(new_phenos) != sum(metadata.phenotype):
                swap_dict = {f:np.random.random()>.5 for f in np.unique(metadata.family)}
                new_phenos = [abs(phenotype-swap_dict[family]) for family, phenotype in zip(metadata.family, metadata.phenotype)]
        elif datatype=='obesity':
            swap_dict = {f:np.random.choice(metadata.phenotype) for f in np.unique(metadata.family)}
            new_phenos = [swap_dict[family] for family in metadata.family]
            while sum(new_phenos) != sum(metadata.phenotype):
                swap_dict = {f:np.random.choice(metadata.phenotype) for f in np.unique(metadata.family)}
                new_phenos = [swap_dict[family] for family in metadata.family]  
        phenos_table[n_iter,:] = new_phenos
        phenos_table[n_iter,:]
    return phenos_table
print('permuting asd data...')
permuted_phenos = Generate_Random_Table(autism_metadata, n_iters=100000, datatype='autism')
np.save(PROJECT_DIR + 'intermediate_files/permutation_test/autism_phenos_permuted.npy', permuted_phenos)

print('permuting obesity data...')
permuted_phenos = Generate_Random_Table(obesity_metadata, n_iters=100000, datatype='obesity')
np.save(PROJECT_DIR + 'intermediate_files/permutation_test/obesity_phenos_permuted.npy', permuted_phenos)


