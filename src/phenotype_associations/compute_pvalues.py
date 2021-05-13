import numpy as np
import pandas as pd
import sys
import scipy.sparse
from tqdm import tqdm
import scipy.stats

cpus = multiprocessing.cpu_count()
BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'

dataset = sys.argv[1]
biomarker = sys.argv[2]
step_size=100000

person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/person_biomarker_%s_%s.npz' % (biomarker, dataset))
sample_data = pd.read_table(BIOMARKER_DIR + 'data/%s/sample_metadata.tsv' % (dataset))
pvalues = np.zeros(np.shape(person_biomarker)[1]) + np.nan
print('person x biomarker shape: ', np.shape(person_biomarker))
for chunk_start in tqdm(np.arange(0,np.shape(person_biomarker)[1], step=step_size)):
    chunk_end = min(chunk_start+step_size, np.shape(person_biomarker)[1])
    person_biomarker_current = person_biomarker[:,chunk_start:chunk_end].toarray()
    affected = person_biomarker_current[sample_data[sample_data.phenotype].index,:]
    unaffected = person_biomarker_current[sample_data[sample_data.phenotype==False].index,:]
    if 'autism' in dataset:
        pvalues[chunk_start:chunk_end] = [scipy.stats.wilcoxon(
            affected[:,i],
            unaffected[:,i]).pvalue for i in range(np.shape(person_biomarker_current)[1])]
    elif 'obesity' in dataset:
        pvalues[chunk_start:chunk_end] = [scipy.stats.mannwhitneyu(
            affected[:,i],
            unaffected[:,i]).pvalue for i in range(np.shape(person_biomarker_current)[1])]
np.save(BIOMARKER_DIR + 'results/phenotype_associations/pvals_%s_%s.npz' % (dataset, biomarker))
print('success -- ', sum(pvalues<(.05/len(pvalues))), 'biomarkers passed bonferonni correction')


