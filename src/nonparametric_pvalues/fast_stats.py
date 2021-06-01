import numpy as np
from scipy.stats import *


def countGreaterThanEqual(pvals_data, pvals_permute):
    idx_data = 0
    idx_permute = 0
    L = len(pvals_permute)
    pvals_count = np.zeros(len(pvals_data)) 
    while (idx_permute<len(pvals_permute)) & (idx_data<len(pvals_data)):
        if pvals_permute[idx_permute]>=pvals_data[idx_data]:
            pvals_count[idx_data] = L-idx_permute
            idx_data = idx_data + 1
        else:
            idx_permute = idx_permute + 1
    return pvals_count

def mannwhitneyu_fast(x, y):
    alternative='two-sided'
    n1 = len(x)
    n2 = len(y)
    u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(x, axis=0)  # calc U for x
    u2 = n1*n2 - u1  # remainder is U for y
    T = tiecorrect(np.concatenate([x,y]))
    if T == 0:
        return 0
    sd = np.sqrt(T * n1 * n2 * (n1+n2+1) / 12.0)
    meanrank = n1*n2/2.0 + 0.5
    bigu = max(u1, u2)
    z = (1*(u1>u2) - 1*(u2>=u1))*(bigu - meanrank) / sd
    return z
    
def wilcoxon_fast(diffs_gt_zero, diffs_rank):
    count = sum(abs(diffs_gt_zero))
    #r_plus = np.sum(x*)
    #r_minus = np.sum(y)
    W = sum(diffs_rank*diffs_gt_zero)
    se = np.sqrt(count*(count+1)*(2*count+1)/6)
    z = W/se
    return z