import pandas as pd
import os
import argparse
from scipy.stats import norm
import json
from multiprocessing import Pool
from functools import partial
from graph import *
import mle
import data
from rpy2.robjects import packages


huge = packages.importr("huge")

p = argparse.ArgumentParser(description="Connect The Dots - Find the most connected subgraph")
p.add_argument("--experimental", help="Experimental dataset file name",
               default='')
p.add_argument("--control", help="Control dataset file name", default='')  # data/example_argininemia/control.csv
p.add_argument("--adj_matrix", help="CSV with adjacency matrix", default='')  # data/example_argininemia/adj.csv
p.add_argument("--s_module",
               help="Comma-separated list or path to CSV of graph G nodes to consider when searching for the most "
                    "connected subgraph.")
p.add_argument("--include_not_in_s",
               help="Include the nodes not appearing in S (encoded with a zero in the optimal bitstring) "
                    "in the most connected subgraph. These are excluded by default.", action="count")
p.add_argument("--kmx", help="Number of highly perturbed nodes to consider. Ignored if S module is given.", default=15,
               type=int)
p.add_argument("--present_in_perc_for_s",
               help="Percentage of patients having metabolite for selection of S module. Ignored if S module is given.",
               default=0.5, type=float)
p.add_argument("--output_name", help="Name of the output JSON file.")
p.add_argument("--out_graph_name", help="Name of the output graph adjecancy CSV file.")
p.add_argument("--glasso_criterion",
               help="Graph-ical Lasso prediction of the graph criterion. 'stars' is default, 'ebic' is faster.",
               default='stars')
p.add_argument("-v", "--verbose", action="count", help="Enable verbose logging.")

if __name__ == '__main__':

    argv = p.parse_args()

    if argv.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Read input dataframe with experimental (positive, disease) samples
    if os.path.exists(argv.experimental):
        experimental_df = pd.read_csv(argv.experimental)
        #   experimental_df[] <- lapply(experimental_df, as.character)  # TODO remove?
        if os.path.exists(argv.control):
            control_data = pd.read_csv(argv.control, index_col=0)
        else:
            control_data = None

        target_patients = list(experimental_df.columns)

        # TODO: Do we need surrogates?
        # Add surrogate disease and surrogate reference profiles based on 1 standard
        # deviation around profiles from real patients to improve rank of matrix when
        # learning Gaussian Markov Random Field network on data. While surrogate
        # profiles are not required, they tend to learn less complex networks
        # (i.e., networks with less edges) and in faster time.

        experimental_df = data.surrogate_profiles(data=experimental_df, ref_data=control_data, std=1)
    else:
        target_patients = []

    # Read input graph (adjacency matrix)
    if os.path.exists(argv.adj_matrix):
        adj_df = pd.read_csv(argv.adj_matrix, index_col=0)  # keep as DataFrame for now
        adj_df.index = adj_df.columns

    else:
        experimental = huge.huge(experimental_df.to_numpy().T, method='glasso')
        # This will take several minutes. For a faster option, you can use the
        # "ebic" criterion instead of "stars", but we recommend "stars".
        logging.debug('Starting huge.select.')
        experimental_select = huge.huge_select(experimental, criterion=argv.glasso_criterion)
        adj_df = [values for (name, values) in experimental_select.items() if name == 'opt.icov'][0]
        np.fill_diagonal(adj_df, 0)
        adj_df = pd.DataFrame(adj_df, index=experimental_df.index, columns=experimental_df.index)
        if argv.out_graph_name:
            adj_df.to_csv(argv.out_graph_name, index=False)

    # The Encoding Process
    G = {}
    for v in adj_df.columns:
        G[v] = 0.0

    # Choose node subset
    kmx = argv.kmx  # Maximum subset size to inspect
    S_set = []

    if os.path.exists(argv.s_module):
        s_module_df = pd.read_csv(argv.s_module)
        S_perturbed_nodes = [str(node) for node in s_module_df.iloc[:, -1]]
    else:
        S_perturbed_nodes = [node.strip() for node in argv.s_module.split(',')]

    logging.debug('Selected perturbed nodes, S = {}'.format(S_perturbed_nodes))

    # Check if all nodes from the s_module are in graph
    for node in S_perturbed_nodes:
        if node not in G:
            logging.debug('Node "{}" not in graph. Exiting program.'.format(node))
            exit(1)

    # Walk through all the nodes in S module
    logging.debug('Get the single-node encoding node ranks starting from each node.')

    pool = Pool() # TODO: add some limit to the number of processes
    ranks_collection = pool.map(partial(single_node_get_node_ranks, G=G, p1=1.0, threshold_diff=0.01, adj_mat=adj_df,
                                        S=S_perturbed_nodes, num_misses=np.log2(len(G))), S_perturbed_nodes)
    pool.close()
    pool.join()

    # Dictionary ranks contains encodings for each node in S_perturbed_nodes
    ranks = {}
    for tup in ranks_collection:
        ranks[tup[1]] = tup[0]

    # Convert to bitstring
    # Get the bitstring associated with the disease module's metabolites
    pt_bs_by_k = mle.get_pt_bs_by_k(S=S_perturbed_nodes, ranks=ranks)

    # Get encoding length of minimum length code word.
    # experimental_df is dataframe with diseases (and surrogates)
    # and z-values for each metabolite
    # TODO: If graph and disease module are given do we still need experimental_df?

    if os.path.exists(argv.experimental):
        data_mx_pvals = experimental_df[target_patients].apply(lambda x: 2 * norm.sf(abs(x)))
        # p-value is area under curve of normal distribution on the right of the
        # specified z-score. sf generates normal distribution with mean=0, std=1,
        # which is exactly what z-scores are
    else:
        data_mx_pvals = np.zeros((1, 1))

    # If we have here specific Patient ID the function will calculate varPvalue
    try:
        pt_id = data_mx_pvals.columns[0]
    except AttributeError:
        pt_id = None

    # Returns a subset of nodes that are highly connected
    res = mle.get_encoding_length(bs=pt_bs_by_k, G=G, pvals=data_mx_pvals.T, pt_id=pt_id)
    ind_mx = np.where(res['d.score'] == res['d.score'].max())
    highest_dscore_paths = res.iloc[ind_mx]

    # Locate encoding (F) with best d-score
    # Tiebreaker 1: If several results have the same d-score take one with longest BS
    ind_F = np.where(highest_dscore_paths['optimalBS'].str.len() ==
                     highest_dscore_paths['optimalBS'].str.len().max())
    ind_F = highest_dscore_paths.iloc[ind_F]

    # Tiebreaker 2: Take the one with the largest subsetSize
    index_highest = np.where(ind_F['subsetSize'] == ind_F['subsetSize'].max())
    ind_F = ind_F.iloc[index_highest]
    ind_F = ind_F.index.values[0]
    F_info = res.iloc[ind_F]

    # You can interpret the probability assigned to this metabolite set by
    # comparing it to a null encoding algorithm, which uses fixed-length codes
    # for all metabolites in the set. The "d.score" is the difference in bitlength
    # Significance theorem, we can estimate the upper bounds on a p-value by 2^-d.score.

    p_value_F = 2.0 ** (-1 * F_info['d.score'])
    # All metabolites in the bitstring
    logging.debug(f'All metabolites in the bitstring: {[d[0] for d in pt_bs_by_k[ind_F]]}')

    # Just the F metabolites that are in S_arg that were were "found"
    keep_nodes = [1]
    if argv.include_not_in_s:
        keep_nodes = [0, 1]

    Fs = [d[0] for d in pt_bs_by_k[ind_F] if d[1] in keep_nodes]

    logging.debug('Set of highly-connected perturbed metabolites F = {} with p-value = {}'.format(Fs, p_value_F))

    kmcm_probability = 2.0 ** (-1 * len(F_info['optimalBS']))
    optimal_bitstring = F_info['optimalBS']

    out_dict = {
        "S_perturbed_nodes": S_perturbed_nodes,
        "F_most_connected_nodes": Fs,
        "p_value": p_value_F,
        "kmcm_probability": kmcm_probability,
        "optimal_bitstring": optimal_bitstring,
        "number_of_nodes_in_G": len(G)
    }

    if not argv.output_name:
        if argv.experimental == '':
            outfname = os.path.basename(argv.adj_matrix).replace('csv', 'json')
        else:
            outfname = os.path.basename(argv.experimental).replace('csv', 'json')
    else:
        outfname = argv.output_name

    with open(outfname, 'w') as f:
        json.dump(out_dict, f, indent=4)

    with open(f'ranks_{outfname}', 'w') as f:
        json.dump(ranks, f, indent=4)
