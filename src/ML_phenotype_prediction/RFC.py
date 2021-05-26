import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glmnet
import sys
import sklearn.model_selection
from sklearn.ensemble import RandomForestClassifier 
import sklearn.tree
import scipy.sparse
import multiprocessing
cpus = multiprocessing.cpu_count()

dataset = sys.argv[1]
biomarker = sys.argv[2]

BIOMARKER_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/'

n_splits=5
n_iters=20
fpr_interp = np.linspace(0,1,100)
tpr_interps = []
aucs = []
person_biomarker = scipy.sparse.load_npz(BIOMARKER_DIR + 'results/generate_biomarkers/sample_vs_biomarker_%s_%s.npz' % (biomarker, dataset))
np.random.seed(42)

for split_idx in range(n_iters):

    kf =  sklearn.model_selection.KFold(n_splits=n_splits, random_state=split_idx)
    idx = [(train_idx, test_idx) for train_idx, test_idx in kf.split([i for i in range(np.shape(person_biomarker)[0])])]
    train_idxs = [i[0] for i in idx]
    test_idxs = [i[1] for i in idx]
    sample_data = pd.read_table(BIOMARKER_DIR + 'data/%s/sample_metadata.tsv' % (dataset))
    phenotype = sample_data.phenotype

    for i in range(n_splits):
        model = RandomForestClassifier(n_jobs=cpus, random_state=i, n_estimators=100) #glmnet.ElasticNet(fit_intercept=False, standardize=False,n_jobs=cpus, alpha=1, random_state=42) #(cv=n_cv, n_jobs=min(cpus,n_cv), selection='random')
        X_train = person_biomarker[train_idxs[i],:]
        y_train = sample_data.phenotype[train_idxs[i]]
        X_test = person_biomarker[test_idxs[i],:]
        y_test = sample_data.phenotype[test_idxs[i]]
        model.fit(X_train.todense(), 1*y_train)
        y_pred = [i[1] for i in model.predict_proba(X_test)]
        fpr, tpr, _ = sklearn.metrics.roc_curve(y_test, y_pred)
        interp_func = scipy.interpolate.interp1d(fpr, tpr)
        tpr_interp = interp_func(fpr_interp)
        tpr_interps = tpr_interps + [list(tpr_interp)]
        auc = sklearn.metrics.roc_auc_score(y_test, y_pred)
        aucs = aucs + [auc]
print(np.mean(aucs))

tpr_df = pd.DataFrame(tpr_interps)
tpr_df.columns = fpr_interp
#tpr_df.columns = ['fpr', 'tpr']
tpr_df.to_csv(BIOMARKER_DIR + 'results/ML_phenotype_prediction/fpr_vs_tpr_%s_%s.tsv' % (biomarker, dataset), sep='\t')
np.savetxt(BIOMARKER_DIR + 'results/ML_phenotype_prediction/auc_%s_%s.txt' % (biomarker, dataset), aucs)