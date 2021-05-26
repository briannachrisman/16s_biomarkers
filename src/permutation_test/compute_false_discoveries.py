import glob
import numpy as np
import pandas as pd
import sys
import scipy.sparse
from tqdm import tqdm
from collections import Counter
import scipy.stats
from matplotlib import pyplot as plt
import sklearn.linear_model
import sklearn.model_selection
import sklearn.metrics
import multiprocessing

cpus = multiprocessing.cpu_count()
BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'

dataset = sys.argv[1]
biomarker = sys.argv[2]

pvals_data = sorted(np.load(BIOMARKER_DIR + 'results/pvalues/pvals_%s_%s.npy' % (dataset, biomarker)))
permuted_phenos = np.load(BIOMARKER_DIR + 'intermediate_files/permutation_test/%s_phenos_permuted.npy' % dataset).astype(bool)
def countEleLessThanOrEqual(arr1, arr2):
    counts = np.zeros(len(arr1))
    for i in range(len(arr1)):
        count = 0
        for j in range(len(arr2)):
            if (arr2[j] < arr1[i]):
                count+= 1
        counts[i] = count
    return counts

pvals_count = np.zeros(len(pvals_data)).astype(int)
step_size=100000
n_iters=10 #len(permuted_phenos)
person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))
print('person x biomarker shape: ', np.shape(person_biomarker))
for chunk_start in tqdm(np.arange(0,np.shape(person_biomarker)[1], step=step_size)):
    chunk_end = min(chunk_start+step_size, np.shape(person_biomarker)[1])
    person_biomarker_current = person_biomarker[:,chunk_start:chunk_end].toarray()
    
    for i in range(n_iters):
        affected = person_biomarker_current[permuted_phenos[i],:]
        unaffected = person_biomarker_current[~permuted_phenos[i],:]
        if 'autism' in dataset:
            pvals_permute = [scipy.stats.wilcoxon(
                affected[:,i],
                unaffected[:,i]).pvalue for i in range(np.shape(person_biomarker_current)[1])]
        elif 'obesity' in dataset:
            pvals_permute = [scipy.stats.mannwhitneyu(
                affected[:,i],
                unaffected[:,i]).pvalue for i in range(np.shape(person_biomarker_current)[1])]
        pvals_count = pvals_count + countEleLessThanOrEqual(pvals_data, pvals_permute)
fdr = [p/(i-1)/n_iters for i,p in enumerate(pvals_count)]
np.save(BIOMARKER_DIR + 'results/permutation_test/fdr_%s_%s.npy' % (dataset, biomarker))
