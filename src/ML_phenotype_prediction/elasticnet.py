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
biomarker = 'sbb2'
dataset = 'autism'
person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))
np.random.seed(42)
n_splits=5
kf =  sklearn.model_selection.KFold(n_splits=n_splits)
idx = [(train_idx, test_idx) for train_idx, test_idx in kf.split([i for i in range(np.shape(person_biomarker)[0])])]
train_idxs = [i[0] for i in idx]
test_idxs = [i[1] for i in idx]
sample_data = pd.read_table(BIOMARKER_DIR + 'data/%s/sample_metadata.tsv' % (dataset))
phenotype = sample_data.phenotype

n_cv=5
y_preds = []
y_tests = []
for i in range(n_splits):
    print(i)
    model = sklearn.linear_model.ElasticNetCV(cv=n_cv, n_jobs=-1, selection='random')
    X_train = person_biomarker[train_idxs[i],:]
    y_train = sample_data.phenotype[train_idxs[i]]
    X_test = person_biomarker[test_idxs[i],:]
    y_test = sample_data.phenotype[test_idxs[i]]
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    y_preds = y_preds + list(y_pred)
    y_tests = y_tests + list(y_test)
    print(sklearn.metrics.roc_auc_score(y_test, y_pred))
np.save(BIOMARKER_DIR + 'results/ML_phenotype_prediction/y_pred_y_true_%s_%s.npz' % (biomarker, dataset), np.array([y_preds, y_tests]))
fpr, tpr, _ = sklearn.metrics.roc_curve(y_tests, y_preds)
aoc = sklearn.metrics.roc_auc_score(y_tests, y_preds)
plt.plot(fpr, tpr)
plt.plot(fpr, fpr, 'r-')
print(aoc)
np.save(BIOMARKER_DIR + 'results/ML_phenotype_prediction/y_pred_y_true_%s_%s.npz' % (biomarker, dataset), np.array([y_preds, y_tests]))