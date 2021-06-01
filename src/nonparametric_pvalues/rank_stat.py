import glob
import numpy as np
import pandas as pd
import sys
import scipy.sparse
from tqdm import tqdm
from collections import Counter
from scipy.stats import rankdata
import scipy.stats.mstats as mstats

sys.path.append('/home/groups/dpwall/briannac/sequence_based_biomarkers/src/nonparametric_pvalues')
from fast_stats import *


BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'

dataset = sys.argv[1]
biomarker = sys.argv[2]

step_size=100000
person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))
sample_data = pd.read_table(BIOMARKER_DIR + 'data/%s/sample_metadata.tsv' % (dataset))
pvalues = np.zeros(np.shape(person_biomarker)[1]) + np.nan


print('person x biomarker shape: ', np.shape(person_biomarker))
for chunk_start in tqdm(np.arange(0,np.shape(person_biomarker)[1], step=step_size)):
    chunk_end = min(chunk_start+step_size, np.shape(person_biomarker)[1])
    person_biomarker_current = person_biomarker[:,chunk_start:chunk_end].toarray()
    if 'autism' in dataset:
        d = person_biomarker_current[::2]-person_biomarker_current[1::2]
        diffs_gt_zeros = 1*(d>0) - 1*(d<0)
        d[d==0] = np.nan
        diffs_rank = mstats.rankdata(abs(np.ma.masked_invalid(d)), axis=0)
        diffs_rank = diffs_rank.transpose()
        diffs_gt_zeros = diffs_gt_zeros.transpose()
        pvalues[chunk_start:chunk_end] = [wilcoxon_fast(d, rank) for d, rank in zip(diffs_gt_zeros, diffs_rank)]
    if 'obesity' in dataset:
        ranks = rankdata(person_biomarker_current, axis=0)
        affected = ranks[np.where(sample_data.phenotype)[0],:]
        unaffected = ranks[np.where(sample_data.phenotype==False)[0],:]
        pvalues[chunk_start:chunk_end] = [mannwhitneyu_fast(
                affected[:,i],
                unaffected[:,i]) for i in range(np.shape(person_biomarker_current)[1])]
np.save(BIOMARKER_DIR + 'results/rank_stats/rank_stat_%s_%s.npy' % (dataset, biomarker), pvalues)
print('success -- min stat', min(pvalues), ' max stat: ', max(pvalues))
