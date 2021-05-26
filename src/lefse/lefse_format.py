import numpy as np
import pandas as pd
import scipy.sparse
BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'

for biomarker_type in ['micropheno4', 'micropheno6', 'micropheno8', 'taxa', 'asv', 'otu90', 'otu95', 'otu97', 'otu99', 'sbb1', 'sbb2', 'sbb3']:
    for dataset in ['obesity', 'autism']:
        print(biomarker_type, dataset)
        sample_vs_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker_type, dataset))
        with open(BIOMARKER_DIR + 'results/generate_biomarkers/biomarker_names_%s_%s.txt' % (biomarker_type, dataset)) as f:
            biomarker_names = [l.replace('\n', '').replace('\t', '.') for l in f.readlines()]
        sample_metadata = pd.read_table((BIOMARKER_DIR + '/results/generate_biomarkers/').replace('results/generate_biomarkers', 'data') + '%s/sample_metadata.tsv' % dataset, index_col=0)
        sample_vs_biomarker_dataframe = pd.DataFrame(sample_vs_biomarker.todense())
        idx = np.argsort(sample_vs_biomarker_dataframe.std().values/sample_vs_biomarker_dataframe.mean().values)[::-1][:1000]
        sample_vs_biomarker_dataframe = sample_vs_biomarker_dataframe[sample_vs_biomarker_dataframe.columns[idx]]
        sample_vs_biomarker_dataframe.columns = [biomarker_names[i] for i in idx]
        sample_vs_biomarker_dataframe.index = sample_metadata.index
        #sample_vs_biomarker_dataframe.insert(0,  'subclass', sample_metadata.subclass)
        sample_vs_biomarker_dataframe.insert(0, 'subject', sample_metadata.subclass)
        sample_vs_biomarker_dataframe.insert(0, 'phenotype', sample_metadata.phenotype==1.0)
        sample_vs_biomarker_dataframe.to_csv(
            BIOMARKER_DIR + 'intermediate_files/lefse/sample_vs_biomarker_lefse_%s_%s.tsv' % (biomarker_type, dataset), sep='\t', index=None)