import numpy as np
import pandas as pd
import scipy.sparse
BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/results/generate_biomarkers/'

for biomarker_type in [#'taxa', 'otu90', 'otu95', 'otu97', 'otu99', 'asv', 'sbb1', 'sbb2', 
                       'sbb3']:
    for dataset in ['autism', 'obesity']:
        print(biomarker_type, dataset)
        sample_vs_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'sample_vs_biomarker_%s_%s.npz' % (biomarker_type, dataset))
        print('loading biomarker names')
        with open(BIOMARKER_DIR + 'biomarker_names_%s_%s.txt' % (biomarker_type, dataset)) as f:
            biomarker_names = [l.replace('\n', '').replace('\t', '.') for l in f.readlines()]
        print('setting up metadata columns...')
        sample_metadata = pd.read_table(BIOMARKER_DIR.replace('results/generate_biomarkers', 'data') + '%s/sample_metadata.tsv' % dataset, index_col=0)
        sample_vs_biomarker_dataframe = pd.DataFrame(sample_vs_biomarker.todense())
        sample_vs_biomarker_dataframe.columns = biomarker_names
        sample_vs_biomarker_dataframe.index = sample_metadata.index
        #sample_vs_biomarker_dataframe.insert(0, 'subject', sample_metadata.subject)
        #sample_vs_biomarker_dataframe.insert(0,  'subclass', sample_metadata.subclass)
        sample_vs_biomarker_dataframe.insert(0, 'phenotype', sample_metadata.phenotype==1.0)
        print('saving...')
        sample_vs_biomarker_dataframe.to_csv(
            BIOMARKER_DIR + 'sample_vs_biomarker_lefse_%s_%s.tsv' % (biomarker_type, dataset), sep='\t')