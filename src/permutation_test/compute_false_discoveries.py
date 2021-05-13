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
START = int(sys.argv[3])
STOP = int(sys.argv[4])
DATA_PVALS_FILE = sys.argv[5]
FD_COUNT_FILE = sys.argv[6]

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

data_pvalues = np.loadtxt(DATA_PVALS_FILE)[START:STOP]

p_values = [[
    -np.log10(wilcoxon(biomarkers[permuted_pheno, i], biomarkers[~permuted_pheno, i]).pvalue) for permuted_pheno in permuted_phenos
] for data_i, i in tqdm.tqdm(enumerate(range(START, min(STOP,np.shape(biomarkers)[1]))))]
fd_counts = [sum([sum(p>=data_p) for p in p_values]) for data_p in sorted(data_pvalues)]
with open(FD_COUNT_FILE, 'w') as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(fd_counts)