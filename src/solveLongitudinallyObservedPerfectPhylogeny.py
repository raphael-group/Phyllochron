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
from itertools import combinations
from scipy.stats import betabinom
from scipy.special import gammaln, logsumexp
from collections import defaultdict

class solveLongitudinallyObservedPerfectPhylogeny():

    def __init__(self, mutation_tree, timepoints = None, df_character_matrix = None, df_total_readcounts = None, df_variant_readcounts = None, snp_list = [],
                 fp = None, fn = None, ado_precision = 15, z = None, threshold = 1, threads = 1, timelimit = 1800, verbose = True, prefix=None, run_pp = False, brute_force = False, use_COMPASS_likelihood = False, compass_assignment=None):
        
        self.mutation_tree = mutation_tree
        
        if prefix is not None:
            self.prefix = prefix

        self.df_character_matrix = df_character_matrix
        self.timepoints = timepoints
        self.ntimepoints = int(max(self.timepoints) + 1)
        self.timepoint_cell_map = defaultdict(list)
        
        self.sample_threshold = np.bincount(timepoints).tolist()
        
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
        if z is None:
            self.z = mutation_tree.shape[0]
        else:
            self.z = z

        self.threshold = threshold
        self.run_pp = run_pp
        self.brute_force = brute_force
        self.use_COMPASS_likelihood = use_COMPASS_likelihood
        
        self.compass_assignment = None
        if compass_assignment is not None:
            self.compass_assignment = compass_assignment
        
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

                    if self.use_COMPASS_likelihood:
                        coeff = solveLongitudinallyObservedPerfectPhylogeny.compute_SNV_loglikelihoods(1, 1, total_reads - variant_reads,  variant_reads, cell_idx, 0.001, 0.001, self.fp + self.fn, 1)[0][0] - solveLongitudinallyObservedPerfectPhylogeny.compute_SNV_loglikelihoods(2, 0, total_reads - variant_reads,  variant_reads, cell_idx, 0.001, 0.001, self.fp + self.fn, 1)[0][0]
                        coeff_mat[cell_idx, mut_idx] = coeff
                    else:
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
        ntimepoints = int(self.ntimepoints)
        nmutations = int(self.nmutations)
        ncells = int(self.ncells)
        nclones = self.mutation_tree.shape[0]

        if self.compass_assignment is not None:
            nclones = self.compass_assignment.shape[1]

               
        mutation_trees = []
        if self.brute_force == True:
            vertices = range(nclones - 1)
            topologies = [nx.DiGraph(edges) for edges in combinations(combinations(vertices, 2), nclones - 2) if nx.is_tree(nx.DiGraph(edges))]
            topologies = [tree for tree in topologies if all(tree.out_degree(n) <= 2 for n in tree.nodes)]
            for top in topologies:
                profiles = np.zeros((nclones, nmutations))
                for n in top.nodes():
                   predecessors = set(top.predecessors(n))  # Use a set for faster lookup
                   profiles[n+1] =  [1 if i in predecessors else 0 for i in range(nmutations)]
                   profiles[n+1][n] = 1
                mutation_trees.append(profiles)            
            
        else:
            mutation_trees.append(self.mutation_tree)
        print(len(mutation_trees))
        print(ntimepoints, nmutations, ncells)
        max_likelihood = -np.inf
        best_solb = None
        best_mutation_tree = None
        prev_solution = None
        for mutation_tree in mutation_trees:
        
            L = np.zeros((nclones, ncells))
            if self.coeff_mat is not None:
                for r in range(nclones):
                    for cell in range(ncells):
                        if self.compass_assignment is not None:
                            L[r][cell] = np.log(max(1e-50, self.compass_assignment[cell,r]))
                        else:
                            L[r][cell] = np.dot(self.coeff_mat[cell], mutation_tree[r])
            else:
                for cell in range(ncells):
                    for r in range(nclones):
                        zero_elem = np.where(self.df_character_matrix.values[cell] == 0, 1, 0)
                        one_elem = np.where(self.df_character_matrix.values[cell] == 1, 1, 0)
                        L[r][cell] = self.fnweight * np.dot(zero_elem, mutation_tree[r]) + self.fpweight * np.dot(one_elem, mutation_tree[r])

            pred_dict = defaultdict(list)
            dir_pred = defaultdict(list)
            for r in range(nclones):
                for r2 in range(nclones):
                    if np.all(mutation_tree[r2][self.mutation_tree[r] == 1] == 1):
                        pred_dict[r2].append(r)
                        if np.linalg.norm(mutation_tree[r2] - mutation_tree[r], ord=1) <= 1.5:
                            dir_pred[r2].append(r)
       
            

            

            model = gp.Model('solveLongitudinallyObservedPerfectPhylogeny')
            model.setParam('TimeLimit', self.worklimit)

            b = model.addVars(ncells, nclones, vtype=gp.GRB.BINARY, name='b')
            tau = model.addVars(nclones, ntimepoints, vtype=gp.GRB.BINARY, name='tau')
            
            tau_pres = model.addVars(nclones, vtype=gp.GRB.BINARY, name='tau_pres')


            ctau = model.addVars(nmutations, ntimepoints, vtype=gp.GRB.BINARY, name='ctau')
            tau_pred  = model.addVars(nclones, ntimepoints, vtype=gp.GRB.BINARY, name='tau_pred')
            
            # Initialize Trivially Longitudinally Observed Assignment

            if prev_solution is not None:
                for i in range(ncells):
                    b[i,0].Start = prev_solution["b"][i,0]
                    for j in range(1,nclones):
                        b[i,j].Start = prev_solution["b"][i,j]
                for t in range(ntimepoints):
                    tau[0,t].Start = prev_solution["tau"][0,t]
                    tau_pred[0,t].Start = prev_solution["tau_pred"][0,t]
                    for j in range(1,nclones):
                        tau[j,t].Start = prev_solution["tau"][j,t]
                        tau_pred[j,t].Start = prev_solution["tau_pred"][j,t]
            else:
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

            for t in range(ntimepoints):
                for j in range(nclones):
                    c_t_sum = gp.LinExpr()
                    for cell in self.timepoint_cell_map[t]:
                        c_t_sum += b[cell,j]
                    model.addConstr(c_t_sum >= int(self.threshold * self.sample_threshold[t]) * tau[j,t])
                    #model.addConstr(c_t_sum >= self.threshold * tau[j,t])

            pres_sum = gp.LinExpr()
            for j in range(nclones):
                pres_sum += tau_pres[j]
            model.addConstr(pres_sum == int(self.z))

            for j in range(nclones):
                for t in range(ntimepoints):
                    for i in self.timepoint_cell_map[t]:
                        model.addConstr(tau[j,t] >= b[i,j])
                        model.addConstr(tau_pres[j] >= b[i,j])

            
            for j in range(nclones):
                pres_bound = gp.LinExpr()
                for t in range(ntimepoints):
                    u_bound = gp.LinExpr()
                    for i in self.timepoint_cell_map[t]:
                        u_bound += b[i,j]
                        pres_bound += b[i,j]
                    model.addConstr(tau[j,t] <= u_bound)
                model.addConstr(tau_pres[j] <= pres_bound)

            if self.run_pp == False:
                
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
            soltau = np.rint(np.reshape(model.getAttr('x', tau).values(), (nclones, ntimepoints)))
            soltaupred = np.rint(np.reshape(model.getAttr('x', tau_pred).values(), (nclones, ntimepoints)))

            prev_solution = {"b": solb, "tau": soltau, "tau_pred": soltaupred}

            self.likelihood = model.ObjVal
            if self.likelihood > max_likelihood:
                best_solb = solb
                best_mutation_tree = mutation_tree
                max_likelihood - self.likelihood

        sol = np.zeros((ncells, nmutations))
        for i in range(ncells):
            for j in range(nclones):
                if best_solb[i,j] > 0:
                    sol[i] = best_mutation_tree[j]

        if self.df_character_matrix is not None:
            df_solb = pd.DataFrame(sol, index = self.df_character_matrix.index, columns = self.df_character_matrix.columns, dtype=int)
        else:
            df_solb = pd.DataFrame(sol, index = self.df_total_readcounts.index, columns = self.df_total_readcounts.columns, dtype=int)

        self.solB = df_solb
        self.solT_mut, self.solT_cell = solveLongitudinallyObservedPerfectPhylogeny.generate_perfect_phylogeny(df_solb, self.prefix, self.timepoints)
        

        print(self.solT_cell)
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
    def generate_perfect_phylogeny(df_binary, prefix, timepoints):
        
        s = []
        t = []

        clone_prev = defaultdict()
        
        
        solT_mut = nx.DiGraph()
        solT_mut.add_node('root')

        solT_cell = nx.DiGraph()
        solT_cell.add_node('root')

        df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    
        
        mapped_mutations = {}
        idx = 1
        for col in sorted(df_binary.columns):
            mapped_mutations[col] = idx
            idx += 1
        mapped_mutations["root"] = 0


        curr_mapping = 0
        for cell_id, row in df_binary.iterrows():
            if cell_id == 'root':
                continue

            curr_node = 'root'
            for idx, column in enumerate(df_binary.columns[row.values == 1]):
                if column not in mapped_mutations:
                    mapped_mutations[column] = curr_mapping
                    curr_mapping += 1
                if curr_node not in mapped_mutations:
                    mapped_mutations[curr_node] = curr_mapping
                    curr_mapping += 1


                if column in solT_mut[curr_node]:
                    curr_node = column

                else:
                    if column in solT_mut.nodes:
                        raise NameError(f'{column} is being repeated')
                    solT_mut.add_edge(curr_node, column)
                    print(curr_node, column)
                    

                    s.append(mapped_mutations[curr_node])
                    t.append(mapped_mutations[column])



                    solT_cell.add_edge(curr_node, column)
                    curr_node = column
            
            solT_cell.add_edge(curr_node, cell_id)
            if mapped_mutations[curr_node] not in clone_prev.keys():
                clone_prev[mapped_mutations[curr_node]] = defaultdict(int)
            clone_prev[mapped_mutations[curr_node]][timepoints[cell_id]] += 1
        
        df_clone_prev = pd.DataFrame.from_dict(clone_prev, orient="index")
        timepoint = []
        clone_id = []
        clone_prev = []
        for i in range(len(mapped_mutations.keys())):
            for j in df_clone_prev.columns:
                timepoint.append(j)
                if i not in df_clone_prev.index:
                    clone_prev.append(0)
                else:
                    clone_prev.append(df_clone_prev.at[i,j])

                clone_id.append(i)

        print(mapped_mutations)
        df_edges = pd.DataFrame({"source": s, "target": t}).drop_duplicates()
        print(df_edges)
        df_edges.to_csv(f'{prefix}_tree_edges.csv', index=False)
        df_prev = pd.DataFrame({"timepoint": timepoint, "clone_id": clone_id, "clonal_prev": clone_prev}) 

        print(df_prev)
        df_prev.to_csv(f'{prefix}_clone_prev.csv', index=False)

        return solT_mut, solT_cell
	
    @staticmethod
    def compute_SNV_loglikelihoods(c_ref, c_alt, r, v, locus, dropout_rate_ref, dropout_rate_alt, sequencing_error_rate, n_cells):
        def log_n_choose_k(n, k):
            return gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1))
        # If homozygous, the copy number of the only allele is irrelevant for the allelic proportion
        if (c_ref==0):
            c_alt=1
        elif (c_alt==0):
            c_ref==1


        discretized_dropout_rate_ref = round(dropout_rate_ref * 1000)
        discretized_dropout_rate_alt = round(dropout_rate_alt * 1000)
        

        dropout_rate_ref = discretized_dropout_rate_ref / 1000.0
        dropout_rate_alt = discretized_dropout_rate_alt / 1000.0

        likelihood_alleles_cells = [[] for _ in range((c_ref + 1) * (c_alt + 1) - 1)]
        dropoutscores = [0.0] * n_cells

        idx = 0
        for k in range(c_ref):
            for l in range(c_alt + 1):
                if k == 0 and l == 0:
                    continue
                            
                seq_error_rates = [0.0] * n_cells
                eps1 = sequencing_error_rate
                eps2 = sequencing_error_rate
                omega = 50 if k == 0 or l == 0 else 8
                

                f = (1.0 * l / (k + l)) * (1 - eps2) + (1.0 * k / (k + l)) * eps1

                for j in range(n_cells):
                    ref_count = r
                    alt_count = v
                    likelihood_dropout = (gammaln(alt_count + omega * f) + 
                                          gammaln(ref_count + omega * (1 - f)) - 
                                          gammaln(alt_count + ref_count + omega) - 
                                          gammaln(omega * f) - 
                                          gammaln(omega * (1 - f)) + 
                                          gammaln(omega))
                    dropoutscores[j] = likelihood_dropout


                likelihood_alleles_cells[idx] = [0.0] * n_cells
                dropout_prob = (log_n_choose_k(c_ref, k) + log_n_choose_k(c_alt, l) + 
                                (c_ref - k) * np.log(dropout_rate_ref) + 
                                k * np.log(1 - dropout_rate_ref) + 
                                (c_alt - l) * np.log(dropout_rate_alt) + 
                                l * np.log(1 - dropout_rate_alt))
                all_dropout_prob = (dropout_rate_ref ** c_ref) * (dropout_rate_alt ** c_alt)
                dropout_prob -= np.log(1 - all_dropout_prob)

                for j in range(n_cells):
                    likelihood_alleles_cells[idx][j] = dropoutscores[j] + dropout_prob

                idx += 1
                
        if c_alt == 0:
            scores_cells = logsumexp(likelihood_alleles_cells[:-1], axis=0)
        else:
            scores_cells = logsumexp(likelihood_alleles_cells[:-2], axis=0)

        dropoutsref = [0.0] * n_cells
        dropoutsalt = [0.0] * n_cells

        idx = 0
        for k in range(c_ref):
            for l in range(c_alt + 1):
                if k == 0 and l == 0:
                    continue

                for j in range(n_cells):
                    config_prob = np.exp(likelihood_alleles_cells[idx][j] - scores_cells[j])
                    dropoutsref[j] += (c_ref - k) * config_prob
                    dropoutsalt[j] += (c_alt - l) * config_prob

                idx += 1

        return scores_cells, dropoutsref, dropoutsalt
    
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
