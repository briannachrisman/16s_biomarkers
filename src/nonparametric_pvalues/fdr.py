import glob
import numpy as np
import pandas as pd
import sys
import scipy.stats.mstats as mstats
import scipy.sparse
from tqdm import tqdm
from collections import Counter
import time
sys.path.append('/home/groups/dpwall/briannac/sequence_based_biomarkers/src/nonparametric_pvalues')
from fast_stats import *


dataset = sys.argv[1]
biomarker = sys.argv[2]
PERMUTE_START = int(sys.argv[3])
PERMUTE_END = int(sys.argv[4])

BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'


pvals_data = sorted(np.load(BIOMARKER_DIR + 'results/pvalues/pvals_%s_%s.npy' % (dataset, biomarker)))
permuted_phenos = np.load(BIOMARKER_DIR + 'intermediate_files/permutation_test/%s_phenos_permuted.npy' % dataset).astype(int)[PERMUTE_START:PERMUTE_END]
pvals_count = np.zeros(len(pvals_data)).astype(int)
step_size=100000
n_iters=len(permuted_phenos)
person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))

t = time.time()
print('person x biomarker shape: ', np.shape(person_biomarker))
for chunk_start in tqdm(np.arange(0,np.shape(person_biomarker)[1], step=step_size)):
    chunk_end = min(chunk_start+step_size, np.shape(person_biomarker)[1])
    person_biomarker_current = person_biomarker[:,chunk_start:chunk_end].toarray()
    if 'autism' in dataset:
        d = person_biomarker_current[::2]-person_biomarker_current[1::2]
        diffs_gt_zeros = 1*(d>0) - 1*(d<0)
        d[d==0] = np.nan
        diffs_rank = mstats.rankdata(np.ma.masked_invalid(d), axis=1)
        diffs_rank = diffs_rank.transpose()
        diffs_gt_zeros = diffs_gt_zeros.transpose()
#        d = (diffs_gt_zero-diffs_lt_zero).transpose()
    if 'obesity' in dataset:
        ranks = rankdata(person_biomarker_current, axis=0)
    for i in range(n_iters):
        if 'autism' in dataset:
            diffs = np.apply_along_axis(arr=d, axis=0,func1d=lambda x: x*permuted_phenos[i,::2])
            pvals_permute = [wilcoxon_fast(rank, diffs_gt_zero) for diffs_gt_zero, rank in zip(diffs_gt_zeros, diffs_rank)]
        elif 'obesity' in dataset:
            affected = ranks[np.where(permuted_phenos[i])[0],:]
            unaffected = ranks[np.where(np.invert(permuted_phenos[i]))[0],:]
            pvals_permute = sorted([mannwhitneyu_fast(
                affected[:,i],
                unaffected[:,i]) for i in range(np.shape(person_biomarker_current)[1])])
        pvals_count = pvals_count + countLessThanEqual(pvals_data, pvals_permute)
#fdr = [p/(i+.0001)/n_iters for i,p in enumerate(pvals_count)]
np.save(BIOMARKER_DIR + 'intermediate_files/nonparametric_pvalues/%s_%s/false_discovery_count_%s_%s_%i.npy' % (dataset, biomarker, dataset, biomarker, PERMUTE_START), pvals_count)
print('saved at ', BIOMARKER_DIR + 'intermediate_files/nonparametric_pvalues/%s_%s/false_discovery_count_%s_%s_%i.npy' % (dataset, biomarker, dataset, biomarker, PERMUTE_START))