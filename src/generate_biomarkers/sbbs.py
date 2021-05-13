import sys
import numpy as np
from scipy.special import comb
from itertools import combinations
import pandas as pd
import scipy.sparse
from tqdm import tqdm
from collections import Counter

order = int(sys.argv[1])
dataset = sys.argv[2]

DATA_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/data'
RESULTS_DIR = '/home/groups/dpwall/briannac/sequence_based_biomarkers/results/generate_biomarkers'

person_asv_filename = '%s/%s/sample_vs_asv.tsv' % (DATA_DIR, dataset)
asv_variant_filename = RESULTS_DIR + '/asv_vs_biomarker_sbb1_%s.npz' % dataset
variant_list_filename = RESULTS_DIR + '/biomarker_names_sbb1_%s.txt' % dataset

asv_vs_biomarker_file_out = RESULTS_DIR + '/asv_vs_biomarker_sbb%i_%s.npz' % (order, dataset)
biomarker_names_file_out = RESULTS_DIR + '/biomarker_names_sbb%i_%s.txt' % (order, dataset)
person_biomarker_file_out = RESULTS_DIR + '/person_biomarker_sbb%i_%s.npz' % (order, dataset)


asv = []
with open(person_asv_filename, 'r') as f:
    people = next(f).strip().split('\t')
    for line in f:
        asv.append(line.strip().split('\t', maxsplit=1)[0])
        
person_asv = np.matrix(pd.read_table(person_asv_filename, index_col=0).transpose())

with open(variant_list_filename, 'r') as f:
    variants = [l.replace('\n', '') for l in f.readlines()]

asv_variant = scipy.sparse.load_npz(asv_variant_filename)
print('Data loaded: ', 'people', len(people), 'asv', len(asv), 'variants', len(variants))

m, n = asv_variant.shape
cached_combs = np.zeros((n, order+1), dtype=int)
for i in range(n):
    cached_combs[i, :] = [comb(i, x) for x in range(order+1)]

r = int(comb(n, order))
print('Number of possible biomarkers: ', r)

# find sparsity pattern
biomarker_exists = np.zeros((r,), dtype=bool)
print("Determining which biomarkers exist...")
for asv_index in tqdm(range(len(asv))):
    var = asv_variant.getrow(asv_index).indices
    p = var.shape[0]
    num_indices = cached_combs[p, order]
    if order > 1:
        new_biomarkers = np.fromiter(combinations(var, order), dtype=np.dtype(','.join(['i']*order)), count=int(comb(p, order))).view(np.dtype('i')).reshape(-1, order)
    else:
        new_biomarkers = var[:, np.newaxis]
    biomarker_indices = (r-1)*np.ones((num_indices,), dtype=int)

    for j in range(order):
        biomarker_indices -= cached_combs[n-1-new_biomarkers[:, j], order-j]
    biomarker_exists[biomarker_indices] = True
    

print('Number of biomarkers that exist:', np.sum(biomarker_exists))


biomarker_mapping = np.cumsum(biomarker_exists)-1
#np.save('%sbiomarker_exists%d' % (output_dir,order), biomarker_exists)

# now calculate person_biomarker
person_biomarker =scipy.sparse.dok_matrix((len(people), np.sum(biomarker_exists)), dtype=np.float64)
asv_biomarker = scipy.sparse.dok_matrix((len(asv), np.sum(biomarker_exists)), dtype=bool)

print("Computing person_vs_biomarker and asv_vs_biomarker...")
for asv_index in tqdm(range(len(asv))):
    #print(asv_index, end=' ')
    var = asv_variant.getrow(asv_index).indices
    p = var.shape[0]
    num_indices = cached_combs[p, order]
    new_biomarkers = np.fromiter(combinations(var, order), dtype=np.dtype(','.join(['i']*order)), count=int(comb(p, order))).view(np.dtype('i')).reshape(-1, order)
    biomarker_indices = (r-1)*np.ones((num_indices,), dtype=int)
    for j in range(order):
        biomarker_indices -= cached_combs[n-1-new_biomarkers[:, j], order-j]
    person_indices = np.where(person_asv[:, asv_index])[0]
    person_biomarker[np.ix_(person_indices, biomarker_mapping[biomarker_indices])] += np.outer(person_asv[person_indices, asv_index], np.ones((num_indices,)))
    asv_biomarker[asv_index,biomarker_mapping[biomarker_indices]] = True
    
def sparse_unique_columns(M):
    M = M.tocsc()
    m, n = M.shape
    if not M.has_sorted_indices:
        M.sort_indices()
    if not M.has_canonical_format:
        M.sum_duplicates()
    sizes = np.diff(M.indptr)
    idx = np.argsort(sizes)
    Ms = M@scipy.sparse.csc_matrix((np.ones((n,)), idx, np.arange(n+1)), (n, n))
    ssizes = np.diff(Ms.indptr)
    ssizes[1:] -= ssizes[:-1]
    grpidx, = np.where(ssizes)
    grpidx = np.concatenate([grpidx, [n]])
    if ssizes[0] == 0:
        counts = [np.array([0, grpidx[0]])]
    else:
        counts = [np.zeros((1,), int)]
    ssizes = ssizes[grpidx[:-1]].cumsum()
    for i, ss in tqdm(enumerate(ssizes)):
        gil, gir = grpidx[i:i+2]
        pl, pr = Ms.indptr[[gil, gir]]
        dv = Ms.data[pl:pr].view(f'V{ss*Ms.data.dtype.itemsize}')
        iv = Ms.indices[pl:pr].view(f'V{ss*Ms.indices.dtype.itemsize}')
        idxi = np.lexsort((dv, iv))
        dv = dv[idxi]
        iv = iv[idxi]
        chng, = np.where(np.concatenate(
            [[True], (dv[1:] != dv[:-1]) | (iv[1:] != iv[:-1]), [True]]))
        counts.append(np.diff(chng))
        idx[gil:gir] = idx[gil:gir][idxi]
    counts = np.concatenate(counts)
    nu = counts.size - 1
    uniques = M@scipy.sparse.csc_matrix((np.ones((nu,)), idx[counts[:-1].cumsum()],
                                   np.arange(nu + 1)), (n, nu))
    return uniques, idx[counts[:-1].cumsum()], counts[1:]

print('Computing biomarkers in LD....')
_, no_ld_idx, _ = sparse_unique_columns(asv_biomarker)
no_ld_idx = sorted(no_ld_idx)

print("Number of biomarkers after LD filtering: ", len(no_ld_idx))

#np.save('%sperson_variant%d_condensed' %  (output_dir,order), person_biomarker)
biomarkers = combinations(range(len(variants)), order)
print('Computing biomarker names...')

#asv_biomarker = scipy.sparse.csr_matrix(asv_biomarker[:,no_ld_idx])
person_biomarker = scipy.sparse.csr_matrix(person_biomarker[:,no_ld_idx])

is_abundant = ((person_biomarker>0).mean(axis=0)>.1).tolist()[0]
print(sum(is_abundant), ' biomarkers > 10% freq')
asv_biomarker = asv_biomarker[:,np.where(is_abundant)[0]]
person_biomarker = person_biomarker[:,np.where(is_abundant)[0]]

unique_idx_current = 0
current_idx = 0
with open(biomarker_names_file_out, 'w') as f:
    for i,biomarker in tqdm(enumerate(biomarkers)):
        if biomarker_exists[i]:  # If biomarker doesn't even exist, just skip.
            if unique_idx_current==len(no_ld_idx):
                break
            elif current_idx==no_ld_idx[unique_idx_current]:
                if is_abundant[unique_idx_current]: 
                    f.write('\t'.join([variants[x] for x in biomarker]) + '\t' + '\t'.join([str(x) for x in biomarker]) + '\n')
                unique_idx_current += 1
            current_idx += 1
print("Writing to files...")
scipy.sparse.save_npz(asv_vs_biomarker_file_out, asv_biomarker)  # Save sparse matrix of ASV x biomarker.
scipy.sparse.save_npz(person_biomarker_file_out, person_biomarker)  # Save sparse matrix of ASV x biomarker.


