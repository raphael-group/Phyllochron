import argparse
import numpy as np
import sys
import pandas as pd




def permute_timepoints(args):

    input_file = args.i
    seed = args.s
    np.random.seed(seed)
    output_file = args.o


    timepoints = pd.read_csv(input_file)
    permuted_timepoints = timepoints.sample(frac=1, random_state=seed).reset_index(drop=True)
    permuted_timepoints.drop("Unnamed: 0", axis=1, inplace=True)
    permuted_timepoints.index.name = ""
    permuted_timepoints.to_csv(output_file)
    return





if __name__ == "__main__":
	

  parser = argparse.ArgumentParser()
  parser.add_argument('-i', type=str, help='input timepoint csv file')
  parser.add_argument('-s', type=int, help='random seed')
  parser.add_argument('-o', type=str, help='output permuted csv file')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  permute_timepoints(args)

