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
        df.to_csv(f'{prefix}gt_imputation_{t}_{error_rate}_{sd}.csv')

        num_to_flip = int(error_rate * character_matrix.size)

        indices_to_flip = np.random.choice(character_matrix.size, size=num_to_flip, replace=False)

        binary_array_flat = character_matrix.flatten()
        binary_array_flat[indices_to_flip] = 1 - binary_array_flat[indices_to_flip]
        character_matrix = binary_array_flat.reshape(character_matrix.shape)

        permuted_character_matrix = character_matrix.copy()


        f = open(f'{prefix}sphyr_{t}_{error_rate}_{sd}.txt', "w")
        f.write(f'{int(t * ncells)}\n')
        f.write(f'{nmutations}\n')
        for r in range(character_matrix.shape[0]):
            f.write("\t".join(map(str,character_matrix[r].flatten())))
            f.write("\n")
        
        
        char_mat = character_matrix[:,:]
        np.savetxt(f'{prefix}scite_{t}_{error_rate}_{sd}.csv', char_mat.T, delimiter=' ', fmt='%d')
        
        timepoints = []
        for tp in range(t):
            for i in range(ncells):
                timepoints.append(tp)
        df = pd.DataFrame(data=timepoints, index=list(range(t * ncells)), columns=['timepoints'])
        df.to_csv(f'{prefix}phyllochron_{t}_{error_rate}_{sd}_timepoints.csv')

        df = pd.DataFrame(data=character_matrix, index=list(range(t * ncells)), columns=mutation_names)
        df.to_csv(f'{prefix}phyllochron_{t}_{error_rate}_{sd}_character_matrix.csv')

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', type=int, help='number of timepoint samples to simulate')
  parser.add_argument('--sd', type=int, help='seed instance for reproducibility')
  parser.add_argument('--ncells', type=int, help='number of cells sequenced per sample')
  parser.add_argument('--error-rate', type=float, help='error rate to simulate under')
  parser.add_argument('--profiles', type=str, help='filepath of the clone profiles csv file')
  parser.add_argument('-o', type=str, help='filepath of output files')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  main(args)
