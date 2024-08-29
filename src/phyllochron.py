#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed August 28 2024

@author: Akhil Jakatdar
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
import networkx as nx

from solveLongitudinallyObservedPerfectPhylogeny import solveLongitudinallyObservedPerfectPhylogeny

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            #if child.startswith('s'):
            subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"

def main(args):
    if args.i is not None:
        df_character_matrix = pd.read_csv(f'{args.i}', index_col = 0)
    if args.r is not None:
        df_total_readcounts = pd.read_csv(f'{args.r}', index_col = 0)
    if args.v is not None:
        df_variant_readcounts = pd.read_csv(f'{args.v}', index_col = 0)
    if args.mutation_tree is not None:
        mutation_tree = pd.read_csv(f'{args.mutation_tree}').values
    if args.t is not None:
        timepoints = pd.read_csv(f'{args.t}', index_col = 0)['timepoints'].values



    fp = args.a
    fn = args.b
    ado = args.ado
    run_pp = args.run_pp

    if args.r is not None:
        solver = solveLongitudinallyObservedPerfectPhylogeny(mutation_tree, timepoints=timepoints, df_total_readcounts=df_total_readcounts,
                                       df_variant_readcounts=df_variant_readcounts, fp=fp, fn=fn,
                                       ado_precision = ado, z=args.z, prefix=args.o, run_pp = run_pp)
    else:
        solver = solveLongitudinallyObservedPerfectPhylogeny(mutation_tree, timepoints=timepoints, df_character_matrix=df_character_matrix,
                                        fp=fp, fn=fn, ado_precision = ado, z=args.z, prefix=args.o, run_pp = run_pp)
        
    solver.solveML_LA()
    prefix = args.o
    solver.writeSolution(f'{prefix}_B.csv')
    solver.writeDOT(f'{prefix}_tree.dot')
    solver.writeDOT(f'{prefix}_tree_without_cells.dot', withcells=False)
    
    if solver.solT_cell is not None:
        with open(f'{prefix}_tree.newick', 'w') as out:
                out.write(tree_to_newick(solver.solT_cell) + ';')
        nx.write_graphml(solver.solT_cell, f'{prefix}_tree.graphml')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='csv file with mutation matrixd')
    parser.add_argument('-r', type=str, help='csv file with total read count matrix')
    parser.add_argument('-v', type=str, help='csv file with variant read count matrix')
    parser.add_argument('-t', type=str, help='csv file with timepoints associated with cells (in the same order that the cells appear in other input files)')
    parser.add_argument('--mutation-tree', type=str, help='csv file with mutation tree encoded by clone profile matrix')
    parser.add_argument('-z', type=float, help='fractional threshold representing the minimum proportion of cells at a sample assigned to a clone for that clone to be present in that sample')
    parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.001)
    parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.001)
    parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=15)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    parser.add_argument('--time-limit', type=int, help='time limit in seconds [1800]', default = 1800)
    parser.add_argument('--run-pp', type=bool, help='run Phyllochron without longitudinal constraints?', default = False)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    if args.mutation_tree is None:
        raise Exception("please provide mutation tree file!")
    if args.t is None and args.run_pp == False:
        raise Exception("please provide timepoint file!")
    if args.i is None and (args.r is None and args.v is None):
        raise Exception("please provide either the binarized character matrix, or the total and variant readcount matrices!")
    
    main(args)
