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
from time import time
import sys

PROJECT_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'


BIOMARKER_FILE =  sys.argv[1]
AUTISM_OR_OBESITY = sys.argv[2]
DATA_PVALS_FILE = sys.argv[3]

print('Loading files...')
biomarkers = np.load(BIOMARKER_FILE)
biomarkers = biomarkers[:,((biomarkers>0).mean(axis=0)>=.1)] # Only include biomarkers in >= 10% of population.

if AUTISM_OR_OBESITY=='obesity':
    person_vs_taxa = pd.read_table(PROJECT_DIR + 'data/obese_lean_twins/person_vs_taxa.tsv')
    metadata =  pd.read_table(PROJECT_DIR + 'data/obese_lean_twins/metadata.tsv', index_col=0)
    metadata = metadata.loc[person_vs_taxa.index]
    metadata['phenotype'] = metadata.obesitycat=='Obese'
    permuted_phenos = np.load(PROJECT_DIR + 'intermediate_files/permutation_test/obesity_phenos_permuted.npy').astype(bool)
    phenos = metadata.phenotype
elif AUTISM_OR_OBESITY=='autism':
    person_vs_taxa = pd.read_table(PROJECT_DIR + 'data/yogurt/person_vs_taxa.tsv')
    metadata =  pd.read_table(PROJECT_DIR + 'data/yogurt/metadata.tsv', index_col=0)
    metadata = metadata.loc[person_vs_taxa.index]
    metadata.phenotype = metadata.phenotype=='A'
    phenos = metadata.phenotype
    permuted_phenos = np.load(PROJECT_DIR + 'intermediate_files/permutation_test/autism_phenos_permuted.npy').astype(bool)

else:
    print('ERROR -- need to specify autism or obesity')
data_pvalues = [-np.log10(wilcoxon(biomarkers[phenos, i], biomarkers[~phenos, i]).pvalue) for i in range(np.shape(biomarkers)[1])]
np.save_txt(DATA_PVALS_FILE, data_pvalues)