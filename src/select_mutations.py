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
        division_result[i] = ad[result[i], :] / dp[result[i], :]
    return division_result

def select_mutations(args):
    patient_name = args.name
    ntimepoints = args.t
    common_mutations = []

    for sample in list(range(1,ntimepoints + 1)):
        print(f'raw_data/{patient_name}-00{sample}.loom')
        ds = lp.connect(f'raw_data/{patient_name}-00{sample}.loom')
        gene=ds.ra.amplicon[0]
        genes = []
        if patient_name == 'AML-88':
            genes = ["FLT3", "FLT3_ITD", "NPM1", "TET2", "DNMT3A", "SF3B1"]
        elif patient_name == 'AML-99':
             genes = ["FLT3","CSF3R", "RUNX1", "DNMT3A", "IDH2", "NRAS", "PTPN11", "IDH1", "SRSF2"]
        elif patient_name == 'AML-63':
             genes = ["FLT3", "FLT3-ITD", "KIT", "NPM1", "IDH2"]
        elif patient_name == 'AML-97':
             genes = ["FLT3", "FLT3-ITD", "WT1", "NPM1", "NRAS", "DNMT3A"]
        elif patient_name == 'AML-38':
             genes = ["FLT3", "FLT3-ITD", "KRAS", "IDH1", "PTPN11", "NPM1", "IDH2", "NRAS"]
        
        A = ahocorasick.Automaton()
        for idx, key in enumerate(genes):
            A.add_word(key, (idx, key))

        A.make_automaton()
        if (sample == 1 and patient_name == 'AML-63') or (sample == 3 and patient_name == "AML-38"):
            common_mutations.append(set(['chr' + ':'.join(j.split(':')[1:]) for i, j in enumerate(ds.ra.id) if len(ds.ra.amplicon[i].split('_')) > 0 and any(A.iter(ds.ra.amplicon[i].split('_')[0]))]))
        elif patient_name == "AML-63" or patient_name == "AML-97" or patient_name == "AML-38":
            common_mutations.append(set([j for i, j in enumerate(ds.ra.id) if len(ds.ra.amplicon[i].split('_')) > 0 and any(A.iter(ds.ra.amplicon[i].split('_')[0]))]))
        else:
            common_mutations.append(set([j for i, j in enumerate(ds.ra.id) if len(ds.ra.amplicon[i].split('_')) > 2 and any(A.iter(ds.ra.amplicon[i].split('_')[1]))]))
    

    common_mutations = list(set.intersection(*common_mutations))

    print(len(common_mutations))
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
        if (sample == 1 and patient_name == 'AML-63') or (sample == 3 and patient_name == 'AML-38'):
            converted_mutations = ['chr' + ':'.join(i.split(':')[1:]) for i in ds.ra.id]
            result = [i for i in range(len(ds.ra.id)) if any(substring in converted_mutations[i] for substring in common_mutations)]

        else:
            result = np.array([i for i, val in enumerate(ds_ra_id_array) if any(A.iter(val))])
            #result = [i for i in range(len(ds.ra.id)) if any(substring in ds.ra.id[i] for substring in common_mutations)]
        
        ds_ad = fast_indexing(ds['AD'][:], result)
        ds_dp = fast_indexing(ds['DP'][:], result)
        vafs = ds_ad/ds_dp
        #vafs = fast_divide(ds['AD'][result,:],ds['DP'][result,:])
        np.nan_to_num(vafs, copy=False,nan=0)
        if patient_name == "AML-99":
            row_indices = np.where(np.sum((vafs >= 0.35), axis=1) >= 0.10*ds.shape[1])[0]
            row_indices_neg = np.where(np.sum((vafs >= 0.35), axis=1) <= 0.90*ds.shape[1])[0]
        elif patient_name == 'AML-88':
            row_indices = np.where(np.sum((vafs >= 0.35), axis=1) >= 0.50*ds.shape[1])[0]
            row_indices_neg = np.where(np.sum((vafs >= 0.35), axis=1) <= 0.50*ds.shape[1])[0]
        elif patient_name == 'AML-63':
            row_indices = np.where(np.sum((vafs >= 0.05), axis=1) >= 0.15*ds.shape[1])[0]
            row_indices_neg = np.where(np.sum((vafs >= 0.05), axis=1) <= 1.00*ds.shape[1])[0]
        elif patient_name == 'AML-97':
            row_indices = np.where(np.sum((vafs >= 0.35), axis=1) >= 0.50*ds.shape[1])[0]
            row_indices_neg = np.where(np.sum((vafs >= 0.35), axis=1) <= 0.50*ds.shape[1])[0]
        elif patient_name == 'AML-38':
            row_indices = np.where(np.sum((vafs >= 0.35), axis=1) >= 0.25*ds.shape[1])[0]
            row_indices_neg = np.where(np.sum((vafs >= 0.35), axis=1) <= 0.75*ds.shape[1])[0]



        if (sample == 1 and patient_name == 'AML-63') or (sample == 3 and patient_name == 'AML-38'):
            mutations_of_interest.append(set([converted_mutations[result[i]] for i in row_indices]))
        else:
            mutations_of_interest.append(set([ds.ra.id[result[i]] for i in row_indices]))
        
        if (sample == 1 and patient_name == 'AML-63') or (sample == 3 and patient_name == 'AML-38'):
            mutations_of_interest_neg.append(set([converted_mutations[result[i]] for i in row_indices_neg]))
        else:
            mutations_of_interest_neg.append(set([ds.ra.id[result[i]] for i in row_indices_neg]))
        
    mutations_of_interest = set.union(*mutations_of_interest)
    mutations_of_interest_neg = set.union(*mutations_of_interest_neg)
    mutations_of_interest = set.intersection(*[mutations_of_interest, mutations_of_interest_neg])
    

    for mutation in sorted(mutations_of_interest):
        print(mutation)
    
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
