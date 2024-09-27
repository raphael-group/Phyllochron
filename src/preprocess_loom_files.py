import loompy
import argparse
import sys
import pandas as pd

def preprocess_loom_files(args):
    patient_name = args.name
    timepoints = args.t
    for timepoint in list(range(1,timepoints + 1)):
        
        print(f'raw_data/{patient_name}-00{timepoint}.loom')
        ds = loompy.connect(f'raw_data/{patient_name}-00{timepoint}.loom')
        allid=ds.ra.id
        df_merged = pd.DataFrame()
        merged_matrix = ds.layer[""][:]

        genes = []
        for i in range(merged_matrix.shape[0]):
            if len(str(ds.ra.amplicon[i]).split('_')) > 1:
                genes.append(str(ds.ra.amplicon[i]).split('_')[1])
            else:
                genes.append(str(ds.ra.amplicon[i]))
        
        special_gene_idxs = [i for i in range(merged_matrix.shape[0]) if genes[i] in ["SRSF2", "RUNX1" , "GATA2"]]
        normal_gene_idxs = [i for i in range(merged_matrix.shape[0]) if genes[i] not in ["SRSF2", "RUNX1" , "GATA2"]]
        
        genotype_layer = ds.layer[""][special_gene_idxs,:]
        ad_layer = ds.layer["AD"][special_gene_idxs,:]
        dp_layer = ds.layer["DP"][special_gene_idxs,:]

        raw_vafs = ad_layer/dp_layer
        raw_vafs[(raw_vafs >= 0.10) & (dp_layer > 5)] = 1
        raw_vafs[(raw_vafs < 0.10) & (dp_layer > 5)] = 0
        raw_vafs[dp_layer < 5] = 3

        merged_matrix[special_gene_idxs,:] = raw_vafs
        
        genotype_layer = ds.layer[""][normal_gene_idxs,:]
        ad_layer = ds.layer["AD"][normal_gene_idxs,:]
        dp_layer = ds.layer["DP"][normal_gene_idxs,:]

        
        raw_vafs = ad_layer/dp_layer
        raw_vafs[(raw_vafs >= 0.10) & (dp_layer > 5)] = 1
        raw_vafs[(raw_vafs < 0.10) & (dp_layer > 5)] = 0
        raw_vafs[dp_layer < 5] = 3
        merged_matrix[normal_gene_idxs,:] = raw_vafs

        outname = f"{args.o}-00{timepoint}.allvariants.genotype_modified.txt"

        df_merged = pd.DataFrame(merged_matrix)
        df_merged.index = allid
        df_merged.to_csv(outname, sep="\t", header=False)

        ds.close()

    for sample in list(range(1,timepoints + 1)):
        ds = loompy.connect(f'raw_data/{patient_name}-00{sample}.loom')
        allid=ds.ra.id
        df_merged = pd.DataFrame()
        merged_matrix = ds.layer[""][:]
    
        genes = []
        for i in range(merged_matrix.shape[0]):
            if len(str(ds.ra.amplicon[i]).split('_')) > 1:
                genes.append(str(ds.ra.amplicon[i]).split('_')[1])
            else:
                genes.append(str(ds.ra.amplicon[i]))
        
        outname = f"processed_loom_files/{patient_name}-00{sample}.variant_readcounts.csv"
        df_merged = pd.DataFrame(ds.layer["AD"][:])
        df_merged.index = allid
        df_merged.to_csv(outname, sep="\t", header=False)

        outname = f"processed_loom_files/{patient_name}-00{sample}.total_readcounts.csv"
        df_merged = pd.DataFrame(ds.layer["DP"][:])
        df_merged.index = allid
        df_merged.to_csv(outname, sep="\t", header=False)
        ds.close()


if __name__ == "__main__":
	

  parser = argparse.ArgumentParser()
  parser.add_argument('-t', type=int, help='number of timepoints sampled')
  parser.add_argument('--name', type=str, help='patient name')
  parser.add_argument('-o', type=str, help='filepath prefix of output genotype file')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  preprocess_loom_files(args)
