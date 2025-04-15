import loompy as lp
import numpy as np
import pickle
import pandas as pd
import argparse 
import sys
from numba import njit, prange

@njit
def fast_indexing(array, indices):
    result = np.empty((len(indices), array.shape[1]))
    for i in range(len(indices)):
        result[i] = array[indices[i], :]
    return result



def produce_SPhyR_input(args):
    patient_name = args.name
    ntimepoints = args.t

    with open(f'mutation_selection/{patient_name}_mutations.pkl', 'rb') as f:
        patient_mutations = pickle.load(f)


    print(len(patient_mutations))
    dfs_to_merge = []
    size = []
    for sample in list(range(1,ntimepoints+1)):
        outname = f"processed_loom_files/{patient_name}-00{sample}.allvariants.genotype_modified.txt"
        print(outname)
        df = pd.read_csv(outname, sep='\t',header=None, index_col=0)
        if sample == 1 and patient_name == 'AML-63':
            df.index = ['chr' + ':'.join(f.split(':')[1:]) for f in df.index]

        filtered_df = df[df.index.str.contains('|'.join(patient_mutations))]
        filtered_index = filtered_df.index
        
    
        filtered_values = filtered_df.values
        zygozity_aware_matrix = np.zeros((filtered_values.shape[0], filtered_values.shape[1]))
        
        nfiltered_index = []
        
        for r in range(filtered_values.shape[0]):
            zygozity_aware_matrix[r, :] = filtered_values[r,:]
            nfiltered_index.append(filtered_index[r] + '_1')
        
        df_merged = pd.DataFrame(zygozity_aware_matrix)
        df_merged.index = nfiltered_index
        dfs_to_merge.append(df_merged)
        size.append(df_merged.shape[1])
    merged_df = pd.concat(dfs_to_merge, axis=1)
    merged_df.fillna(-1, inplace=True)
    merged_df.replace(3.0, -1, inplace=True)


    merged_val = merged_df.values.astype(int)
    merged_val = merged_val.T
    f = open(f'method_input/{patient_name}_sphyr.csv', "w")
    f.write(f'{merged_val.shape[0]}\n')
    f.write(f'{merged_val.shape[1]}\n')
    for r in range(merged_val.shape[0]):
        f.write(" ".join(map(str,merged_val[r].flatten())))
        f.write("\n")

    return merged_df, size

def produce_readcount_input(args):
    patient_name = args.name
    ntimepoints = args.t

    with open(f'mutation_selection/{patient_name}_mutations.pkl', 'rb') as f:
        patient_mutations = pickle.load(f)


    
    print('got here')
    
    dfs_to_merge = []
    size = []
    for sample in list(range(1,ntimepoints+1)):
        outname = f"processed_loom_files/{patient_name}-00{sample}.total_readcounts.csv"
        df = pd.read_csv(outname, sep='\t',header=None, index_col=0)
        if sample == 1 and patient_name == 'AML-63':
            df.index = ['chr' + ':'.join(f.split(':')[1:]) for f in df.index]

        filtered_df = df[df.index.str.contains('|'.join(patient_mutations))]
        filtered_index = filtered_df.index
        filtered_values = filtered_df.values
        zygozity_aware_matrix = np.zeros((filtered_values.shape[0], filtered_values.shape[1]))
    
        nfiltered_index = []
        

        for r in range(filtered_values.shape[0]):
            zygozity_aware_matrix[r, :] = filtered_values[r,:]
            nfiltered_index.append(filtered_index[r] + '_1')
        
        df_merged = pd.DataFrame(zygozity_aware_matrix)
        df_merged.index = nfiltered_index
        dfs_to_merge.append(df_merged)
        size.append(df_merged.shape[1])
    merged_df = pd.concat(dfs_to_merge, axis=1).transpose()   # Initialize with the first DataFrame
    merged_df.fillna(-1, inplace=True)
    merged_df.index = range(len(merged_df))

    
    timepoints = []

    data = np.zeros((ntimepoints, len(merged_df.columns)))
    data_total = np.zeros((ntimepoints, len(merged_df.columns)))

    for tp in range(ntimepoints):
        for i in range(size[tp]):
            timepoints.append(tp)
    df = pd.DataFrame(data=timepoints, index=list(range(len(merged_df))), columns=['timepoints'])
    df.to_csv(f'method_input/{patient_name}_timepoints.csv')


    for i, c in enumerate(merged_df.columns):
        for r in merged_df.index:
            if merged_df.loc[r,c] < 5:
                data[timepoints[r], i] += 1
            data_total[timepoints[r], i] += 1


    df_missing = pd.DataFrame(data/data_total, columns=merged_df.columns)
    #selected_mutations = df_missing.loc[:, ((df_missing >= 0.05)).any()].columns
    selected_mutations = df_missing.columns

    all_zero_indices = merged_df[selected_mutations].index[
        (merged_df[selected_mutations] == 0).sum(axis=1) >= (len(selected_mutations))
    ].tolist()
    merged_df[selected_mutations].drop(all_zero_indices, axis=0).reset_index(drop=True).to_csv(f'method_input/{patient_name}_total_readcounts.csv')    
    
    df.drop(all_zero_indices, axis=0).reset_index(drop=True).to_csv(f'method_input/{patient_name}_timepoints.csv')


    dfs_to_merge = []
    size = []
    for sample in list(range(1,ntimepoints+1)):
        outname = f"processed_loom_files/{patient_name}-00{sample}.variant_readcounts.csv"
        df = pd.read_csv(outname, sep='\t',header=None, index_col=0)
        if sample == 1 and patient_name == 'AML-63':
            df.index = ['chr' + ':'.join(f.split(':')[1:]) for f in df.index]
        

        filtered_df = df[df.index.str.contains('|'.join(patient_mutations))]
        filtered_index = filtered_df.index
        
        #if sample == 1:
            #filtered_index = ['chr' + ':'.join(f.split(':')[1:]) for f in filtered_df.index]
    
        filtered_values = filtered_df.values
        zygozity_aware_matrix = np.zeros((filtered_values.shape[0], filtered_values.shape[1]))
    
        nfiltered_index = []
        
        for r in range(filtered_values.shape[0]):
            zygozity_aware_matrix[r, :] = filtered_values[r,:]
            #zygozity_aware_matrix[r, :][zygozity_aware_matrix[2*r, :] == 2] = 0
            #zygozity_aware_matrix[2*r+1, :][zygozity_aware_matrix[2*r+1, :] == 1] = 0
            #zygozity_aware_matrix[r, :][zygozity_aware_matrix[r, :] == 2] = 1
            nfiltered_index.append(filtered_index[r] + '_1')
            #nfiltered_index.append(filtered_index[r] + '_2')
        
        df_merged = pd.DataFrame(zygozity_aware_matrix)
        df_merged.index = nfiltered_index
        #df_merged = df_merged.loc[:, ~(df_merged == 3).any()]
        dfs_to_merge.append(df_merged)
        size.append(df_merged.shape[1])
        #print(df)
    merged_df_var = pd.concat(dfs_to_merge, axis=1).transpose()   # Initialize with the first DataFrame
    merged_df_var.fillna(-1, inplace=True)
    merged_df_var.index = range(len(merged_df_var))
    merged_df_var[selected_mutations].drop(all_zero_indices, axis=0).reset_index(drop=True).to_csv(f'method_input/{patient_name}_variant_readcounts.csv')    


    num_cells = merged_df.drop(all_zero_indices, axis=0).reset_index(drop=True).shape[0]
    with open(f'method_input/{patient_name}_COMPASS_variants.csv', "w") as vcf:
        vcf.write("CHR,REF,ALT,REGION,NAME,FREQ," + ",".join([str(i) for i in range(num_cells)]) + "\n")            
        for i, mutation in enumerate(selected_mutations):
            print(mutation)
            chrom = mutation.split(':')[0]
            pos = mutation
            vid = '.'

            vs = list(merged_df_var.drop(all_zero_indices, axis=0).reset_index(drop=True)[mutation])
            ts = list(merged_df.drop(all_zero_indices, axis=0).reset_index(drop=True)[mutation])
            ref = 'A'
            alt = 'T'
            qual = '.'
            filter_ = '.'
            info = '.'


            vcf.write(f"{chrom},{ref},{alt},{pos},{pos},0," + ",".join([f'{int(tot-var)}:{int(var)}' for var,tot in zip(vs, ts)])+"\n")


    with open(f'method_input/{patient_name}_data_size.pkl', 'wb') as f:
        pickle.dump(size, f)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', type=int, help='number of timepoints sampled')
  parser.add_argument('--name', type=str, help='patient name')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  produce_SPhyR_input(args)
  produce_readcount_input(args)
