from collections import defaultdict
import numpy as np
import scipy.stats
from scipy.stats import betabinom
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import seaborn as sns
import argparse
import sys
import math


def generate_perfect_phylogeny_matrix(graph, nmutations):
    index_map = {}
    root = [node for node, indegree in graph.in_degree() if indegree == 0][0] 
    normal_profile = np.zeros(nmutations) 
    profiles = [normal_profile]
    def rec_generation(curr_node, curr_profile, profiles):
        index_map[curr_node] = len(profiles)
        profiles.append(curr_profile)
        for child in list(graph.successors (curr_node)): 
            new_profile = curr_profile.copy() 
            new_profile[int(child)] = 1
            rec_generation (child, new_profile, profiles)

    f_profile = normal_profile.copy()
    f_profile[int(root)] = 1
    rec_generation (root, f_profile, profiles)
    matrix = np.array(profiles)
    return matrix, index_map

def main(args):
    t = args.t
    error_rate = args.error_rate
    read_depth = args.rd
    ado = args.ado
    ncells = args.ncells
    nmutations = args.nmutations
    prefix = args.o
    sd = args.sd
    np.random.seed(args.sd)

    mutation_names = [f'm{i}' for i in range(nmutations)]
    #random_tree = nx.random_tree(nmutations)
    #directed_tree = nx.bfs_tree(random_tree, root)
    random_tree = nx.DiGraph()
    random_tree.add_node(0)  # Add root node
    # Grow the tree by preferentially attaching new nodes to existing nodes
    for i in range(1, nmutations):
        # Pick a random existing node with a bias towards nodes with low indegree (favoring depth)
        # Here, you could bias towards nodes with a smaller degree for deeper trees
        candidates = [node for node in random_tree.nodes if random_tree.out_degree(node) < 2]
        if candidates:
            parent = np.random.choice(candidates)
            random_tree.add_edge(parent, i)
    
    def remove_nodes_keep_connected_digraph(tree, p):
        """
        Randomly remove nodes from a directed tree with probability p while keeping it connected.
        """
        G = tree.copy()
        for node in list(G.nodes):
            # Skip root and leaf nodes to avoid breaking connectivity
            if G.in_degree(node) == 1 and G.out_degree(node) > 0 and np.random.random() < p:
                parent = next(G.predecessors(node))
                children = list(G.successors(node))
                
                # Connect each child to the node's parent
                for child in children:
                    G.add_edge(parent, child)
                    
                # Remove the node
                G.remove_node(node)
        
        return G
    

    print(len(random_tree.nodes()))
    print(random_tree.edges())
    profiles,index_map = generate_perfect_phylogeny_matrix(random_tree, nmutations)

    random_tree = remove_nodes_keep_connected_digraph(random_tree, 0.5)

    # Continuous Tumor Growth Simulations
    max_depth = max([sum(profiles[i]) for i in range(profiles.shape[0])])
    depth_increments = math.floor(max_depth/t)
    continuous_growth_profile_map = defaultdict(list)
    for tp in range(t - 1):
        for i in range(profiles.shape[0]):
            if sum(profiles[i]) < (tp+1) * depth_increments:
                continuous_growth_profile_map[tp].append(i)
    continuous_growth_profile_map[t - 1] = [i for i in range(profiles.shape[0])]
    print(continuous_growth_profile_map)

    # Sliding Window Lifespan  Simulations
    sliding_window_profile_map = defaultdict(list)
    for tp in range(t):
        for i in range(profiles.shape[0]):
            if tp < t - 1 and sum(profiles[i]) < (tp+1) * depth_increments and sum(profiles[i]) >= tp * depth_increments:
                sliding_window_profile_map[tp].append(i)
            elif tp == t - 1 and sum(profiles[i]) >= tp * depth_increments:
                sliding_window_profile_map[tp].append(i)

    print(sliding_window_profile_map)

    # Random trivial graph reduction
    trivial_reduction_profile_map = defaultdict(list)
    
    random_lifespan_tree = nx.DiGraph()
    for n in random_tree.nodes():
        for tp in range(t):
            random_lifespan_tree.add_node(f"{n}_{tp}")
        for tp in range(1,t):
            random_lifespan_tree.add_edge(f"{n}_{tp-1}", f"{n}_{tp}")

    for (u,v) in random_tree.edges():
        for tp in range(t):
            random_lifespan_tree.add_edge(f"{u}_{tp}", f"{v}_{tp}")
        for tp in range(t - 1):
            random_lifespan_tree.add_edge(f"{u}_{tp}", f"{v}_{tp+1}")

    # Randomly delete edges w/ probability p s.t. resulting graph is LC
    def delete_nodes_with_connectivity(G, p):
        # Make a copy of the graph to avoid modifying the original
        G_copy = G.copy()
        
        # Keep trying until we get a single connected component after node deletions
        while True:
            # Generate a list of nodes to delete with probability p
            nodes_to_delete = [node for node in G_copy.nodes if np.random.random() < p]
            
            # Remove the selected nodes
            G_temp = G_copy.copy()
            G_temp.remove_nodes_from(nodes_to_delete)
            
            isValid = True
            
            for n in range(nmutations):
                for t1 in range(t - 2):
                    for t2 in range(t1 + 1, t - 1):
                        for t3 in range(t2 + 1, t):
                            if f"{n}_{t1}" in G_temp.nodes() and f"{n}_{t2}" not in G_temp.nodes() and f"{n}_{t3}" in G_temp.nodes():
                                isValid = False
            

            def is_reachable_from_single_root(G):
                # Find nodes with in-degree 0
                roots = [node for node in G.nodes if G.in_degree(node) == 0]

                # There must be exactly one node with in-degree 0
                if len(roots) != 1:
                    return False

                # Get the single root node
                root = roots[0]

                # Check if all nodes are reachable from the root node
                reachable_nodes = nx.descendants(G, root)  # All nodes reachable from root
                reachable_nodes.add(root)  # Include the root itself

                return len(reachable_nodes) == len(G)

            # Check if the resulting graph is strongly connected
            if isValid and len(nodes_to_delete) > 0 and is_reachable_from_single_root(G_temp):
                return G_temp  # Return the modified graph if it's still strongly connected
            


    print(len(random_lifespan_tree.nodes()))
    random_lifespan_tree = delete_nodes_with_connectivity(random_lifespan_tree, p=0.50)

    print(len(random_lifespan_tree.nodes()))
    
    for clone in random_tree.nodes():
        for tp in range(t):
            if f"{clone}_{tp}" in random_lifespan_tree.nodes():
                trivial_reduction_profile_map[tp].append(index_map[clone])


    print(trivial_reduction_profile_map)
    character_matrix = np.zeros((t * ncells, nmutations))

    for tp in range(t):
        n_clones = len(trivial_reduction_profile_map[tp]) + 1
        proportion_rare = 0.05
        if n_clones != 1:
            
            # Randomly select one clone to be rare
            low_prob_idx = np.random.choice(n_clones)

            # Sample random probabilities for all clones
            prob_clones = np.random.rand(n_clones)

            # Set chosen low-probability clone to exactly 0.1
            prob_clones[low_prob_idx] = 0.0  # temporarily set to zero

            # Apply minimum threshold to remaining clones
            for idx in range(n_clones):
                if idx != low_prob_idx:
                    prob_clones[idx] = max(prob_clones[idx], proportion_rare)

            # Normalize remaining clones to sum to (1 - 0.1)
            prob_clones /= prob_clones.sum()
            prob_clones *= (1.0 - proportion_rare)

            # Now create the full probability vector
            full_prob_clones = prob_clones.copy()
            full_prob_clones[low_prob_idx] = proportion_rare
            print(full_prob_clones)
            # Calculate how many cells to assign to rare clone
            n_low_prob_cells = int(proportion_rare * ncells)
            n_remaining_cells = ncells - n_low_prob_cells

            # Assign rare clone
            assignments = [low_prob_idx] * n_low_prob_cells
        else:
            full_prob_clones = [1.0]
            n_remaining_cells = ncells
            assignments = []
        # Assign remaining clones
        remaining_assignments = np.random.choice(
            np.arange(n_clones),
            size=n_remaining_cells,
            p=full_prob_clones
        )

        assignments.extend(remaining_assignments)
        np.random.shuffle(assignments)

        # Assign profiles
        for i, clone_idx in enumerate(assignments):
            character_matrix[tp * ncells + i] = profiles[clone_idx]
        '''
        prob_clones = np.maximum(np.random.rand(len(trivial_reduction_profile_map[tp]) + 1), 0.05)
        normalized_prob_clones = prob_clones/prob_clones.sum()
        for i in range(ncells):
            #character_matrix[tp * ncells + i] = profiles[np.random.choice(continuous_growth_profile_map[tp])]
            character_matrix[tp * ncells + i] = profiles[np.random.choice([0] + trivial_reduction_profile_map[tp], p=normalized_prob_clones)]
        '''
    og_character_matrix = character_matrix.copy()
    
    df_profiles = pd.DataFrame(data=profiles.astype(int), index=list(range(nmutations + 1)), columns=mutation_names)
    df_profiles.to_csv(f'{prefix}gt_tree_{ncells}_{t}_{error_rate}_{sd}.csv', sep='\t', index=False)

    df = pd.DataFrame(data=og_character_matrix, index=list(range(t * ncells)), columns=mutation_names)
    df.to_csv(f'{prefix}gt_readcount_{ncells}_{t}_{error_rate}_{sd}.csv')
    
    total_rc_matrix = character_matrix.copy()
    variant_rc_matrix = character_matrix.copy()

    for r in range(og_character_matrix.shape[0]):
        for c in range(og_character_matrix.shape[1]):
            total_rc_matrix[r,c] = int(np.random.poisson(read_depth))
            fp_error_rate = 0.001 + (1 - 0.001) * 0.50 * og_character_matrix[r,c]
            alpha = fp_error_rate * ado
            beta = (1-fp_error_rate) * ado
            variant_rc_matrix[r,c] = int(scipy.stats.betabinom.rvs(int(total_rc_matrix[r,c]), alpha, beta))
    
    nmissing = math.floor(error_rate * og_character_matrix.shape[0] * og_character_matrix.shape[1])
    selected_cell_indices = np.random.randint(og_character_matrix.shape[0], size=nmissing)
    

    #selected_character_indices = np.random.randint(og_character_matrix.shape[1], size=nmissing)
    selected_character_indices = np.random.randint(low=3, high=7, size=nmissing)
    total_rc_matrix[selected_cell_indices, selected_character_indices] = 0
    variant_rc_matrix[selected_cell_indices, selected_character_indices] = 0
    df_total = pd.DataFrame(data=total_rc_matrix, index=list(range(t * ncells)), columns=mutation_names)
    df_total.to_csv(f'{prefix}phyllochron_{ncells}_{t}_{error_rate}_{sd}_total_readcounts.csv')

    df_variant = pd.DataFrame(data=variant_rc_matrix, index=list(range(t * ncells)), columns=mutation_names)
    df_variant.to_csv(f'{prefix}phyllochron_{ncells}_{t}_{error_rate}_{sd}_variant_readcounts.csv')
    
    
    with open(f'{prefix}phylovar_{ncells}_{t}_{error_rate}_{sd}.mpileup', "w") as out:
        str_muts = [str(i) for i in range(nmutations)]

        for i, mutation in enumerate(df_variant.columns):
            vs = list(df_variant[mutation])
            ts = list(df_total[mutation])
            
            seq_name = f'chr1'
            pos = i+1
            ref_allele = 'A'
            entries = [seq_name, pos, ref_allele]
            
            for var,tot in zip(vs, ts):
                if tot == 0:
                    tot = 1
                readcov = int(tot)
                readstring = '.'*int(tot-var) + 'T'*int(var)
                qualstring = 'I'*int(tot)
                entries += [readcov, readstring, qualstring]
            out.write("\t".join(map(str, entries)))
            if i < nmutations - 1:
                out.write('\n')
                
    
    with open(f'{prefix}phylovar_{ncells}_{t}_{error_rate}_{sd}.names', "w") as out:
        for i in range(ncells*t):
            out.write(f'{i}\tCT\n')
            
    with open(f'{prefix}compass_{ncells}_{t}_{error_rate}_{sd}_variants.csv', "w") as vcf:
        
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
    

    VAF_mat = variant_rc_matrix/total_rc_matrix
    mutation_mat = ((VAF_mat >= 0.1) & (variant_rc_matrix >= 3)).astype(int)
    mutation_mat[selected_cell_indices, selected_character_indices] = -1
    print(mutation_mat.shape)

    f = open(f'{prefix}sphyr_{ncells}_{t}_{error_rate}_{sd}.txt', "w")
    f.write(f'{int(t * ncells)}\n')
    f.write(f'{nmutations}\n')
    for r in range(mutation_mat.shape[0]):
        f.write("\t".join(map(str,mutation_mat[r].flatten())))
        f.write("\n")
    mutation_mat[selected_cell_indices, selected_character_indices] = 3
    char_mat = mutation_mat[:,:]
    np.savetxt(f'{prefix}scite_{ncells}_{t}_{error_rate}_{sd}.csv', char_mat.T, delimiter=' ', fmt='%d')

    np.savetxt(f'{prefix}sclongtree_{ncells}_{t}_{error_rate}_{sd}_genotype_matrix.csv', char_mat, delimiter='\t', fmt='%d')

    timepoints = []
    for tp in range(t):
        for i in range(ncells):
            timepoints.append(tp)

    with open(f'{prefix}sclongtree_{ncells}_{t}_{error_rate}_{sd}_timepoints.csv', "w") as cell_timepoints:
        for tp in range(t):
            cell_idxs = ";".join([str(idx) for idx, val in enumerate(timepoints) if val == tp])
            np.savetxt(f'{prefix}sclongtree_{ncells}_{t}_{error_rate}_{sd}_{tp}_genotype_matrix.csv', char_mat[[int(idx) for idx, val in enumerate(timepoints) if val == tp]], delimiter='\t', fmt='%d')
            cell_timepoints.write(f"t{tp}\t{cell_idxs}\n")
    
    
    df = pd.DataFrame(data=timepoints, index=list(range(t * ncells)), columns=['timepoints'])
    df.to_csv(f'{prefix}phyllochron_{ncells}_{t}_{error_rate}_{sd}_timepoints.csv')


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', type=int, help='number of timepoint samples to simulate')
  parser.add_argument('-rd', type=int, help='average sequencing read depth')
  parser.add_argument('--ado', type=int, help='ADO precision parameter')
  parser.add_argument('--sd', type=int, help='seed instance for reproducibility')
  parser.add_argument('--ncells', type=int, help='number of cells sequenced per sample')
  parser.add_argument('--error-rate', type=float, help='error rate to simulate under')
  parser.add_argument('--profiles', type=str, help='filepath of the clone profiles csv file')
  parser.add_argument('--nmutations', type=int, help='number of somatic mutations sequenced')
  parser.add_argument('-o', type=str, help='filepath of output files')

  args = parser.parse_args(None if sys.argv[1:] else ['-h']) 
  main(args)
