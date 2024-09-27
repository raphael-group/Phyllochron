import loompy as lp
import numpy as np
import pickle
import pandas as pd
import argparse 
import sys
import ahocorasick
from numba import njit, prange
import random

@njit
def fast_indexing(array, indices):
    result = np.empty((len(indices), array.shape[1]))
    for i in range(len(indices)):
        result[i] = array[indices[i], :]
    return result

@njit(parallel=True)
def fast_divide(ad, dp, result):
    division_result = np.empty((len(result), ad.shape[1]))
    for i in prange(len(result)):
        print(i)
        division_result[i] = ad[result[i], :] / dp[result[i], :]
    return division_result

def select_mutations(args):
    patient_name = args.name
    ntimepoints = args.t
    common_mutations = []

    for sample in list(range(1,ntimepoints + 1)):
        print(f'raw_data/{patient_name}-00{sample}.loom')
        ds = lp.connect(f'raw_data/{patient_name}-00{sample}.loom')
        print(ds.shape)
        if sample == 1 and patient_name == 'AML-63':
            common_mutations.append(set(['chr' + ':'.join(i.split(':')[1:]) for i in ds.ra.id]))
        else:
            common_mutations.append(set(ds.ra.id))

    for i in range(len(common_mutations)):
        if 'chr5:170837543:C/CTCTG' in common_mutations[i]:
            print('NPM1', i)
    for i in range(len(common_mutations)):
        if 'chr1:115258744:C/T' in common_mutations[i]:
            print('NRAS', i)
    common_mutations = list(set.intersection(*common_mutations))
    print(common_mutations[0])
    if 'chr5:170837543:C/CTCTG' in common_mutations:
        print('NPM1 found')

    
    print('common_mutations', len(common_mutations))

    mutations_of_interest = []
    mutations_of_interest_neg = []
    A = ahocorasick.Automaton()
    for idx, key in enumerate(common_mutations):
        A.add_word(key, (idx, key))

    A.make_automaton()
    for sample in list(range(1,ntimepoints + 1)):
        
        print(f'raw_data/{patient_name}-00{sample}.loom')
        ds = lp.connect(f'raw_data/{patient_name}-00{sample}.loom')
        ds_ra_id_array = np.array(ds.ra.id)
        if sample == 1 and patient_name == 'AML-63':
            converted_mutations = ['chr' + ':'.join(i.split(':')[1:]) for i in ds.ra.id]
            result = [i for i in range(len(ds.ra.id)) if any(substring in converted_mutations[i] for substring in common_mutations)]

        else:
            result = np.array([i for i, val in enumerate(ds_ra_id_array) if any(A.iter(val))])
            #result = [i for i in range(len(ds.ra.id)) if any(substring in ds.ra.id[i] for substring in common_mutations)]
        
        print('result', len(result))
        ds_ad = fast_indexing(ds['AD'][:], result)
        ds_dp = fast_indexing(ds['DP'][:], result)
        vafs = ds_ad/ds_dp
        #vafs = fast_divide(ds['AD'][result,:],ds['DP'][result,:])
        print('got here')
        np.nan_to_num(vafs, copy=False,nan=0)
        row_indices = np.where(np.sum((vafs >= 0.35), axis=1) >= 0.4*ds.shape[1])[0]

        row_indices_neg = np.where(np.sum((vafs >= 0.35), axis=1) <= 0.5*ds.shape[1])[0]
        

        if sample == 1 and patient_name == 'AML-63':
            mutations_of_interest.append(set([converted_mutations[result[i]] for i in row_indices]))
        else:
            mutations_of_interest.append(set([ds.ra.id[result[i]] for i in row_indices]))

        if sample == 1 and patient_name == 'AML-63':
            mutations_of_interest_neg.append(set([converted_mutations[result[i]] for i in row_indices_neg]))
        else:
            mutations_of_interest_neg.append(set([ds.ra.id[result[i]] for i in row_indices_neg]))
        
    mutations_of_interest = set.union(*mutations_of_interest)
    mutations_of_interest_neg = set.union(*mutations_of_interest_neg)
    mutations_of_interest = set.intersection(*[mutations_of_interest, mutations_of_interest_neg])
    
    random.seed(42)
    if len(mutations_of_interest) > 15:
        mutations_of_interest = random.sample(mutations_of_interest, 15)
    annotations = {i:"_".join(i.split(':')) for i in mutations_of_interest}
    
    with open(f'{args.o}_mutations.pkl', 'wb') as f:
        pickle.dump(mutations_of_interest, f)

    with open(f'{args.o}_annotations.pkl', 'wb') as f:
        pickle.dump(annotations, f)

if __name__ == "__main__":
	

  parser = argparse.ArgumentParser()
  parser.add_argument('-t', type=int, help='number of timepoints sampled')
  parser.add_argument('--name', type=str, help='patient name')
  parser.add_argument('-o', type=str, help='filepath prefix of output selected mutation file')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  select_mutations(args)
