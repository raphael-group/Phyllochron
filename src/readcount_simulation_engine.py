import numpy as np
import scipy.stats
from scipy.stats import betabinom
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import argparse
import sys

def main(args):
    t = args.t
    error_rate = args.error_rate
    read_depth = args.rd
    ado = args.ado
    ncells = args.ncells
    profiles = pd.read_csv(f'{args.profiles}').values
    mutation_names = pd.read_csv(f'{args.profiles}').columns
    nmutations = profiles.shape[1]
    prefix = args.o
    sd = args.sd
    np.random.seed(sd)

    character_matrix = np.zeros((t * ncells, nmutations))

    for tp in range(t):
        for i in range(ncells):
            if tp < t//2:
                j = np.random.randint(3)
                character_matrix[tp * ncells + i] = profiles[j]
            if tp == t//2:
                j = np.random.randint(3)
                if j == 0:
                    character_matrix[tp * ncells + i] = profiles[0]
                else:
                    character_matrix[tp * ncells + i] = profiles[3]

            if tp >= t//2 + 1:
                j = np.random.randint(4,7)
                if j > 5:
                    character_matrix[tp * ncells + i] = profiles[0]
                else:
                    character_matrix[tp * ncells + i] = profiles[j]

    og_character_matrix = character_matrix.copy()
    df = pd.DataFrame(data=og_character_matrix, index=list(range(t * ncells)), columns=mutation_names)
    df.to_csv(f'{prefix}gt_readcount_{t}_{error_rate}_{sd}.csv')

    total_rc_matrix = character_matrix.copy()
    variant_rc_matrix = character_matrix.copy()

    for r in range(og_character_matrix.shape[0]):
        for c in range(og_character_matrix.shape[1]):
            total_rc_matrix[r,c] = int(np.random.poisson(read_depth))
            fp_error_rate = error_rate/2 + (1 - error_rate) * 0.50 * og_character_matrix[r,c]
            alpha = fp_error_rate * ado
            beta = (1-fp_error_rate) * ado
            variant_rc_matrix[r,c] = int(scipy.stats.betabinom.rvs(int(total_rc_matrix[r,c]), alpha, beta))
    
    df_total = pd.DataFrame(data=total_rc_matrix, index=list(range(t * ncells)), columns=mutation_names)
    df_total.to_csv(f'{prefix}phyllochron_{t}_{error_rate}_{sd}_total_readcounts.csv')

    df_variant = pd.DataFrame(data=variant_rc_matrix, index=list(range(t * ncells)), columns=mutation_names)
    df_variant.to_csv(f'{prefix}phyllochron_{t}_{error_rate}_{sd}_variant_readcounts.csv')
    
    
    with open(f'{prefix}phylovar_{t}_{error_rate}_{sd}.mpileup', "w") as out:
        str_muts = [str(i) for i in range(nmutations)]

        for i, mutation in enumerate(df_variant.columns):
            vs = list(df_variant[mutation])
            ts = list(df_total[mutation])
            
            seq_name = f'chr1'
            pos = i+1
            ref_allele = 'A'
            entries = [seq_name, pos, ref_allele]
            
            for var,tot in zip(vs, ts):
                readcov = int(tot)
                readstring = '.'*int(tot-var) + 'T'*int(var)
                qualstring = 'I'*int(tot)
                entries += [readcov, readstring, qualstring]
            out.write("\t".join(map(str, entries)))
            if i < nmutations - 1:
                out.write('\n')
                
    
    with open(f'{prefix}phylovar_{t}_{error_rate}_{sd}.names', "w") as out:
        for i in range(ncells*t):
            out.write(f'{i}\tCT\n')
            
    with open(f'{prefix}compass_{t}_{error_rate}_{sd}_variants.csv', "w") as vcf:
        
        # Write header line
        vcf.write("CHR,REF,ALT,REGION,NAME,FREQ," + ",".join([str(i) for i in range(t*ncells)]) + "\n")            
        for i, mutation in enumerate(df_variant.columns):
            chrom = 'chr1'
            pos = f'region{i+1}'
            vid = '.'
            
            vs = list(df_variant[mutation])
            ts = list(df_total[mutation])
            ref = 'A'
            alt = 'T'
            qual = '.'
            filter_ = '.'
            info = '.'
            
            
            vcf.write(f"{chrom},{ref},{alt},{pos},{pos},0," + ",".join([f'{tot-var}:{var}' for var,tot in zip(vs, ts)])+"\n")
    timepoints = []
    for tp in range(t):
        for i in range(ncells):
            timepoints.append(tp)
    
    df = pd.DataFrame(data=timepoints, index=list(range(t * ncells)), columns=['timepoints'])
    df.to_csv(f'{prefix}phyllochron_{t}_{error_rate}_{sd}_timepoints.csv')


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', type=int, help='number of timepoint samples to simulate')
  parser.add_argument('-rd', type=int, help='average sequencing read depth')
  parser.add_argument('--ado', type=int, help='ADO precision parameter')
  parser.add_argument('--sd', type=int, help='seed instance for reproducibility')
  parser.add_argument('--ncells', type=int, help='number of cells sequenced per sample')
  parser.add_argument('--error-rate', type=float, help='error rate to simulate under')
  parser.add_argument('--profiles', type=str, help='filepath of the clone profiles csv file')
  parser.add_argument('-o', type=str, help='filepath of output files')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  main(args)
