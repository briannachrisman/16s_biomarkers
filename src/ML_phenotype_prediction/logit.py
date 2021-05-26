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
import scipy.interpolate
from glmnet import LogitNet
import seaborn as sns

dataset = sys.argv[1]
biomarker = sys.argv[2]


cpus = multiprocessing.cpu_count()
BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'


n_splits=5
n_cv=5

person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))
np.random.seed(42)

kf =  sklearn.model_selection.KFold(n_splits=n_splits)
idx = [(train_idx, test_idx) for train_idx, test_idx in kf.split([i for i in range(np.shape(person_biomarker)[0])])]
train_idxs = [i[0] for i in idx]
test_idxs = [i[1] for i in idx]
sample_data = pd.read_table(BIOMARKER_DIR + 'data/%s/sample_metadata.tsv' % (dataset))
phenotype = sample_data.phenotype

y_preds = []
y_tests = []
fpr_interp = np.linspace(0,1,100)
tpr_interps = []
aocs = []
for i in range(n_splits):
    print(i)
    model = LogitNet(fit_intercept=False, n_jobs=cpus) #(cv=n_cv, n_jobs=min(cpus,n_cv), selection='random')
    X_train = person_biomarker[train_idxs[i],:]
    y_train = sample_data.phenotype[train_idxs[i]]
    X_test = person_biomarker[test_idxs[i],:]
    y_test = sample_data.phenotype[test_idxs[i]]
    model.fit(X_train.todense(), 1*y_train)
    y_pred = model.predict_proba(X_test)[:,1]
    fpr, tpr, _ = sklearn.metrics.roc_curve(y_test, y_pred)
    interp_func = scipy.interpolate.interp1d(fpr, tpr)
    tpr_interp = interp_func(fpr_interp)
    tpr_interps = tpr_interps + [list(tpr_interp)]
    aoc = sklearn.metrics.roc_auc_score(y_test, y_pred)
    aocs = aocs + [aoc]
tpr_df = pd.DataFrame(tpr_interps)
tpr_df.columns = fpr_interp
tpr_df = tpr_df.melt()
tpr_df.columns = ['fpr', 'tpr']
tpr_df.to_csv(BIOMARKER_DIR + 'results/ML_phenotype_prediction/fpr_vs_tpr_%s_%s.tsv' % (biomarker, dataset), sep='\t')
print(np.mean(aocs))
