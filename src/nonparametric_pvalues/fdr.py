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
SEED_IDX = int(sys.argv[3])

BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'


pvals_data = sorted(np.load(BIOMARKER_DIR + 'results/rank_stats/rank_stat_%s_%s.npy' % (dataset, biomarker)))
permuted_phenos = np.load(BIOMARKER_DIR + 'intermediate_files/permutation_test/%s_phenos_permuted.npy' % dataset).astype(int)
#permuted_phenos[:,::2] = 0
#permuted_phenos[:,1::2] = 1

pvals_count = np.zeros(len(pvals_data)).astype(int)
biomarker_sample_size = 100000
n_biomarker_samples = 10
n_iters=1
person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))
print(np.shape(person_biomarker))
pvals_total = []
t = time.time()
for i_biomarker_sample in range(n_biomarker_samples):
    print(i_biomarker_sample)
    np.random.seed(i_biomarker_sample + 10000*SEED_IDX)
    person_biomarker_current = person_biomarker[:,np.random.choice([i for i in range(np.shape(person_biomarker)[1])], biomarker_sample_size, replace=True)].toarray()
    if 'autism' in dataset:
        d = person_biomarker_current[::2]-person_biomarker_current[1::2]
        diffs_gt_zeros = 1*(d>0) - 1*(d<0)
        d[d==0] = np.nan
        diffs_rank = mstats.rankdata(abs(np.ma.masked_invalid(d)), axis=0)
        diffs_rank = diffs_rank.transpose()
        diffs_gt_zeros = diffs_gt_zeros.transpose()
    if 'obesity' in dataset:
        ranks = rankdata(person_biomarker_current, axis=0)
    for i in range(n_iters): # MAKE RAND_IDX A CHOICE
        rand_idx=np.random.choice([i for i in range(np.shape(person_biomarker_current)[1])])
        if 'autism' in dataset:
            diffs_gt_zero_permute = diffs_gt_zeros*(2*permuted_phenos[rand_idx,::2]-1)
            pvals_permute = [wilcoxon_fast(diffs_gt_zero, rank) for diffs_gt_zero, rank in zip(diffs_gt_zero_permute, diffs_rank)]
        elif 'obesity' in dataset:
            affected = ranks[np.where(permuted_phenos[rand_idx])[0],:]
            unaffected = ranks[np.where(permuted_phenos[rand_idx]==False)[0],:]
            pvals_permute = sorted([mannwhitneyu_fast(
                affected[:,i],
                unaffected[:,i]) for i in range(np.shape(person_biomarker_current)[1])])
        #pvals_count = pvals_count + countGreaterThanEqual(pvals_data, pvals_permute)/len(pvals_permute)/n_iters/n_biomarker_samples]
        pvals_total = pvals_total + list(pvals_permute)
print('max pval permute: ', max(pvals_total))

#fdr = [p/(i+.0001)/n_iters for i,p in enumerate(pvals_count)]
#np.save(BIOMARKER_DIR + 'intermediate_files/nonparametric_pvalues/%s_%s/false_discovery_count_%s_%s_%i.npy' % (dataset, biomarker, dataset, biomarker, SEED_IDX), pvals_count)
np.savetxt(BIOMARKER_DIR + 'intermediate_files/nonparametric_pvalues/%s_%s/permuted_rank_stat_%s_%s_%i.txt' % (dataset, biomarker, dataset, biomarker, SEED_IDX), pvals_total)

print('saved at ', BIOMARKER_DIR + 'intermediate_files/nonparametric_pvalues/%s_%s/false_discovery_count_%s_%s_%i.npy' % (dataset, biomarker, dataset, biomarker, SEED_IDX))
