#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed August 28 2024

@author: Akhil Jakatdar
"""

import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools
from scipy.stats import betabinom
from collections import defaultdict

class solveLongitudinallyObservedPerfectPhylogeny():

    def __init__(self, mutation_tree, timepoints = None, df_character_matrix = None, df_total_readcounts = None, df_variant_readcounts = None, snp_list = [],
                 fp = None, fn = None, ado_precision = 15, z = 0.05, threads = 1, timelimit = 1800, verbose = True, prefix=None, run_pp = False):
        
        self.mutation_tree = mutation_tree
        
        if prefix is not None:
            self.prefix = prefix

        self.df_character_matrix = df_character_matrix
        self.timepoints = timepoints
        self.ntimepoints = int(max(self.timepoints) + 1)
        self.timepoint_cell_map = defaultdict(list)
        
        self.threshold = np.bincount(timepoints).tolist()
        
        for i in range(len(self.timepoints)):
            self.timepoint_cell_map[self.timepoints[i]].append(i)
        
        # read count matrices
        self.df_total_readcounts = df_total_readcounts        
        self.df_variant_readcounts = df_variant_readcounts        

        if df_character_matrix is not None:
            self.mutation_list = list(df_character_matrix.columns)
            self.ncells = len(self.df_character_matrix)
            self.nmutations = len(self.df_character_matrix.columns)
        else:
            self.mutation_list = list(df_total_readcounts.columns)
            self.ncells = len(self.df_total_readcounts)
            self.nmutations = len(self.df_total_readcounts.columns)

        self.fp = fp
        self.fn = fn
        self.z = z
        self.run_pp = run_pp
        
        
        if df_total_readcounts is not None:
            self.cell_list = list(df_total_readcounts.index)
            self.mutation_list = list(df_total_readcounts.columns)

            bb_alpha = fp * ado_precision
            bb_beta = (1 - fp) * ado_precision
            
            coeff_mat = np.zeros((self.ncells, self.nmutations))
            for cell_idx, cell in enumerate(self.cell_list):
                for mut_idx, mutation in enumerate(self.mutation_list):
                    total_reads = df_total_readcounts.loc[cell][mutation]
                    variant_reads = df_variant_readcounts.loc[cell][mutation]
                    if total_reads > 0:
                        coeff = betabinom.logpmf(variant_reads, total_reads, 1, 1) - betabinom.logpmf(variant_reads, total_reads, bb_alpha, bb_beta)
                        coeff_mat[cell_idx, mut_idx] = coeff
            self.coeff_mat = coeff_mat                    
        else:
            self.coeff_mat = None
        
        # gurobi parameters
        self.threads = threads
        self.worklimit = timelimit
        self.verbose = verbose
        self.run_pp = run_pp
        
        self.fpweight = np.log(1 - fp) - np.log(fp)
        self.fnweight = np.log(fn) - np.log(1 - fn)
        
        # solution
        # self.B = np.zeros((self.ncells, self.nmutations))
        self.solB = None
        self.solT_cell = None
        self.solT_mut = None

    def solveML_LA(self):
        '''
        first_self.mutation_tree = np.array([[0,0,0,0,0,0],
            [0,0,0,0,1,0],
            [0,0,0,0,1,1],
            [0,0,1,0,1,1],
            [0,1,1,0,1,1],
            [0,1,1,1,1,1],
            [1,1,1,0,1,1],
            ])

        first_self.mutation_tree = np.array([[0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 1, 0, 1],
            [0, 0, 1, 0, 0, 1, 1, 1],
            [0, 0, 1, 0, 1, 1, 1, 1],
            [0, 1, 1, 0, 1, 1, 1, 1],
            [1, 1, 1, 0, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1]]
           )
        first_self.mutation_tree = np.array([[0,0,0,0,0],
            [0,0,0,0,1],
            [0,0,0,1,1],
            [0,0,1,0,1],
            [0,1,1,0,1],
            [1,0,1,0,1],
            ])


        '''

        ntimepoints = int(self.ntimepoints)
        nmutations = int(self.nmutations)
        ncells = int(self.ncells)
        nclones = self.mutation_tree.shape[0]

        print(ntimepoints, nmutations, ncells)
        L = np.zeros((nclones, ncells))
        if self.coeff_mat is not None:
            for r in range(nclones):
                for cell in range(ncells):
                    L[r][cell] = np.dot(self.coeff_mat[cell], self.mutation_tree[r])
        else:
            for cell in range(ncells):
                for r in range(nclones):
                    zero_elem = np.where(self.df_character_matrix.values[cell] == 0, 1, 0)
                    one_elem = np.where(self.df_character_matrix.values[cell] == 1, 1, 0)
                    L[r][cell] = self.fnweight * np.dot(zero_elem, self.mutation_tree[r]) + self.fpweight * np.dot(one_elem, self.mutation_tree[r])

        pred_dict = defaultdict(list)
        dir_pred = defaultdict(list)
        for r in range(nclones):
            for r2 in range(nclones):
                if np.all(self.mutation_tree[r2][self.mutation_tree[r] == 1] == 1):
                    pred_dict[r2].append(r)
                    if np.linalg.norm(self.mutation_tree[r2] - self.mutation_tree[r], ord=1) <= 1.5:
                        dir_pred[r2].append(r)
   
        

        

        model = gp.Model('solveLongitudinallyObservedPerfectPhylogeny')
        model.setParam('TimeLimit', self.worklimit)

        b = model.addVars(ncells, nclones, vtype=gp.GRB.BINARY, name='b')
        tau = model.addVars(nclones, ntimepoints, vtype=gp.GRB.BINARY, name='tau')
        ctau = model.addVars(nmutations, ntimepoints, vtype=gp.GRB.BINARY, name='ctau')
        tau_pred  = model.addVars(nclones, ntimepoints, vtype=gp.GRB.BINARY, name='tau_pred')
        
        # Initialize Trivially Longitudinally Observed Assignment
        for i in range(ncells):
            b[i,0].Start = 1
            for j in range(1,nclones):
                b[i,j].Start = 0
        for t in range(ntimepoints):
            tau[0,t].Start = 1
            tau_pred[0,t].Start = 1
            for j in range(1,nclones):
                tau[j,t].Start = 0
                tau_pred[j,t].Start = 0

        for i in range(ncells):
            assignment_cons = gp.LinExpr()
            for j in range(nclones):
                assignment_cons += b[i,j]

            model.addConstr(assignment_cons == 1)


        if self.run_pp == False:
            for t in range(ntimepoints):
                for j in range(nclones):
                    c_t_sum = gp.LinExpr()
                    for cell in self.timepoint_cell_map[t]:
                        c_t_sum += b[cell,j]
                    model.addConstr(c_t_sum >= int(self.z * self.threshold[t]) * tau[j,t])
                    # real data
                    #15
                    #30
                    #400
            for j in range(nclones):
                for t in range(ntimepoints):
                    for i in self.timepoint_cell_map[t]:
                        model.addConstr(tau[j,t] >= b[i,j])
            
            for j in range(nclones):
                for t in range(ntimepoints):
                    u_bound = gp.LinExpr()
                    for i in self.timepoint_cell_map[t]:
                        u_bound += b[i,j]

                    model.addConstr(tau[j,t] <= u_bound)

            for t1 in range(ntimepoints - 2):
                for t2 in range(t1 + 1, ntimepoints - 1):
                    for t3 in range(t2 + 1, ntimepoints):
                        for j in range(nclones):
                            model.addConstr(tau[j,t1] + (1 - tau[j,t2]) + tau[j,t3] <= 2)
            
            for t in range(ntimepoints):
                for j in range(nclones):
                    pred_cons = gp.LinExpr()
                    for pred in pred_dict[j]:
                        model.addConstr(tau_pred[j,t] >= tau[pred,t])
                        pred_cons += tau[pred,t]

                    model.addConstr(tau_pred[j,t] <= pred_cons)


            for t1 in range(ntimepoints - 2):
                for t2 in range(t1, ntimepoints - 1):
                    for j in range(nclones):
                        for pred in pred_dict[j]:
                            model.addConstr(tau[j, t1] + (1 - tau[pred, t2]) + tau[pred, t2 + 1] <= 2)
            
            for i in range(ncells):
                t = self.timepoints[i]
                for j in range(nclones):
                    for t1 in range(ntimepoints):
                        if t1 <= t:
                            model.addConstr(tau_pred[j,t1] >= b[i,j])

        obj_sum = gp.LinExpr()

        for i in range(ncells):
            for j in range(nclones):
                obj_sum += L[j, i] * b[i,j]

        
        model.setObjective(obj_sum, gp.GRB.MAXIMIZE)
        
        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.setParam(gp.GRB.Param.Method, 4)
        
        model.setParam(gp.GRB.Param.FeasibilityTol, 1e-6)
        model.setParam(gp.GRB.Param.IntFeasTol, 1e-6)
        model.setParam(gp.GRB.Param.OptimalityTol, 1e-6)
        
        model.optimize()
        solb = np.rint(np.reshape(model.getAttr('x', b).values(), (ncells, nclones)))
        sol = np.zeros((ncells, nmutations))

        self.likelihood = model.getObjective()
        for i in range(ncells):
            for j in range(nclones):
                if solb[i,j] > 0:
                    sol[i] = self.mutation_tree[j]

        if self.df_character_matrix is not None:
            df_solb = pd.DataFrame(sol, index = self.df_character_matrix.index, columns = self.df_character_matrix.columns, dtype=int)
        else:
            df_solb = pd.DataFrame(sol, index = self.df_total_readcounts.index, columns = self.df_total_readcounts.columns, dtype=int)

        self.solB = df_solb
        self.solT_mut, self.solT_cell = solveLongitudinallyObservedPerfectPhylogeny.generate_perfect_phylogeny(df_solb)
        
        sol = df_solb.values
        tsol = np.hstack((sol, self.timepoints.reshape(-1, 1)))
        print(np.unique(tsol, axis=0))
        
        return
        
    def writeSolution(self, fname):
        if self.solB is not None:
            self.solB.to_csv(fname)
    
    @staticmethod
    def expand_multi_state_to_binary(df_multistate):

        ncells = len(df_multistate)
        binarized_mat = None
        binary_col_dict = {}
        for column in df_multistate.columns:
            max_state = df_multistate[column].max()
            for s in range(1, max_state+1):
                state_col = np.zeros((ncells))
                if s == 1:
                    state_col[df_multistate[column] > 0] = 1
                else:
                    state_col[df_multistate[column] == s] = 1

                binary_col_dict[f'{column}_{s}'] = state_col

        df_binary = pd.DataFrame(binary_col_dict, index = df_multistate.index, dtype=int)
        return df_binary    
    
    @staticmethod
    def generate_perfect_phylogeny(df_binary):
        

        solT_mut = nx.DiGraph()
        solT_mut.add_node('root')

        solT_cell = nx.DiGraph()
        solT_cell.add_node('root')

        df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    

        for cell_id, row in df_binary.iterrows():
            if cell_id == 'root':
                continue

            curr_node = 'root'
            for column in df_binary.columns[row.values == 1]:
                if column in solT_mut[curr_node]:
                    curr_node = column
                else:
                    if column in solT_mut.nodes:
                        raise NameError(f'{column} is being repeated')
                    solT_mut.add_edge(curr_node, column)
                    solT_cell.add_edge(curr_node, column)
                    curr_node = column

            solT_cell.add_edge(curr_node, cell_id)   

        return solT_mut, solT_cell

    def writeDOT(self, dot_file, withcells=True):
        if withcells is True:
            writeTree = self.solT_cell
        else:
            writeTree = self.solT_mut
        
        with open(dot_file, 'w') as output:

            output.write(f'digraph N {{\n')
            output.write(f"\toverlap=\"false\"\n")
            output.write(f"\trankdir=\"TB\"\n")

            idx_dict = {}
            idx = 0
            if writeTree is not None:
                for node in writeTree.nodes:
                    idx_dict[node] = idx
                    output.write(f'\t{idx} [label=\"{node}\", style=\"bold\"];\n')
                    idx += 1

                for edge in writeTree.edges:
                    output.write(f"\t{idx_dict[edge[0]]} -> {idx_dict[edge[1]]} [style=\"bold\"];\n")

            output.write(f'}}')
