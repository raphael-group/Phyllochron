import loompy as lp
import numpy as np
import pickle
import pandas as pd
import argparse 
import sys



def produce_SPhyR_input(args):
    patient_name = args.name
    ntimepoints = args.t

    with open(f'mutation_selection/{patient_name}_mutations.pkl', 'rb') as f:
        patient_mutations = pickle.load(f)



    dfs_to_merge = []
    size = []
    for sample in list(range(1,ntimepoints+1)):
        outname = f"/n/fs/ragr-data/users/aj7381/{patient_name}-00{sample}.allvariants.genotype_modified.txt"
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
    merged_df_var.to_csv(f'method_input/{patient_name}_variant_readcounts.csv')    


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
    merged_df.to_csv(f'method_input/{patient_name}_total_readcounts.csv')    

    timepoints = []
    for tp in range(ntimepoints):
        for i in range(size[tp]):
            timepoints.append(tp)
    df = pd.DataFrame(data=timepoints, index=list(range(len(merged_df))), columns=['timepoints'])
    df.to_csv(f'method_input/{patient_name}_timepoints.csv')

    with open(f'method_input/{patient_name}_COMPASS_variants.csv', "w") as vcf:
        vcf.write("CHR,REF,ALT,REGION,NAME,FREQ," + ",".join([str(i) for i in range(merged_df.shape[0])]) + "\n")            
        for i, mutation in enumerate(merged_df.columns):
            print(mutation)
            chrom = mutation.split(':')[0]
            pos = mutation
            vid = '.'

            vs = list(merged_df_var[mutation])
            ts = list(merged_df[mutation])
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
