import argparse
import logging
import os
from collections import Counter
import numpy as np
import pandas as pd
from igraph import Graph
from rpy2.robjects import r, pandas2ri, packages, globalenv, ListVector, StrVector
from scipy.stats import norm
import json

huge = packages.importr("huge")
r.source("./CTD/R/data.surrogateProfiles.r")
r.source("./CTD/R/data.imputeData.r")
r.source("./CTD/R/singleNode.getNodeRanksN.r")
r.source("./CTD/R/graph.diffuseP1.r")
r.source("./CTD/R/graph.connectToExt.r")
r.source("./CTD/R/stat.fishersMethod.r")
r.source("./CTD/R/mle.getEncodingLength.r")
r.source("./CTD/R/mle.getPtBSbyK.r")
r.source("./getNodeNames.R")

pandas2ri.activate()

p = argparse.ArgumentParser(description="Connect The Dots - Find the most connected sub-graph")
p.add_argument("--experimental", help="Experimental dataset file name", default = '')  # data/example_argininemia/experimental.csv
p.add_argument("--control", help="Control dataset file name", default = '')            # data/example_argininemia/control.csv
p.add_argument("--adj_matrix", help="CSV with adjacency matrix", default = '')         # data/example_argininemia/adj.csv
p.add_argument("--s_module", help="Comma-separated list or path to CSV of graph G nodes to consider when searching for the most connected sub-graph.")
p.add_argument("--kmx", help="Number of highly perturbed nodes to consider. Ignored if S module is given.", default=15, type=int)
p.add_argument("--present_in_perc_for_s", help="Percentage of patients having metabolite for selection of S module. Ignored if S module is given.", default=0.5, type=float)
p.add_argument("--output_name", help="Name of the output JSON file.")
p.add_argument("--out_graph_name", help="Name of the output graph adjecancy CSV file.")
p.add_argument("--glasso_criterion", help="Graph-ical Lasso prediction of the graph criterion. 'stars' is default, 'ebic' is faster.", default='stars')
p.add_argument("-v", "--verbose", action="count", help="Enable verbose logging.")
argv = p.parse_args()

if argv.verbose:
    logging.getLogger().setLevel(logging.DEBUG)

# Read input dataframe with experimental (positive, disease) samples
if os.path.exists(argv.experimental):
    experimental_df = pd.read_csv(argv.experimental)
    #   experimental_df[] <- lapply(experimental_df, as.character)  # TODO remove?
    control_data = pd.read_csv(argv.control, index_col=0)
    target_patients = list(experimental_df.columns)

    # Prepare for calling R function
    data_surrogateProfiles = globalenv['data.surrogateProfiles']

    experimental_df = data_surrogateProfiles(data=experimental_df, std=1, ref_data=control_data)
else:
    target_patients = np.nan

# Read input graph (adjacency matrix)
if os.path.exists(argv.adj_matrix):
    adj_df = pd.read_csv(argv.adj_matrix)  # keep as DataFrame for now
    adj_df.index = adj_df.columns
else:
    experimental = huge.huge(experimental_df.to_numpy().T, method='glasso')
    # This will take several minutes. For a faster option, you can use the
    # "ebic" criterion instead of "stars", but we recommend "stars".
    logging.debug('Starting huge.select.')
    experimental_select = huge.huge_select(experimental, criterion=argv.glasso_criterion)
    adj_df = [values for (name, values) in experimental_select.items() if name == 'opt.icov'][0]
    np.fill_diagonal(adj_df, 0)
    adj_df = pd.DataFrame(adj_df, index=experimental_df.index, columns=experimental_df.index)  # keep as DataFrame for now
    if argv.out_graph_name:
        adj_df.to_csv(argv.out_graph_name, index=False)

logging.debug('Convert adjacency matrices to an igraph object.')
igraph = Graph.Weighted_Adjacency(adj_df, mode='max')  # TODO: potential optimization

# The Encoding Process
adj_mat = igraph.get_adjacency(attribute="weight")

G = {}
for v in igraph.vs:
    G[v['name']] = 0

# Choose node subset
kmx = argv.kmx  # Maximum subset size to inspect
S_set = []

if not argv.s_module:
    for pt in target_patients:
        temp = experimental_df.sort_values(by=pt, ascending=False)
        S_patient = list(temp.index)[:kmx]
        S_set += S_patient

    # Created list containing top kmx metabolites for every target user
    occurrences = Counter(S_set)

    # Keep in the S module the metabolites perturbed in at least 50% patients
    S_perturbed_nodes = [node for node in occurrences if occurrences[node] >= len(target_patients) * argv.present_in_perc_for_s]


elif os.path.exists(argv.s_module):
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

# Prepare for calling R function
singleNode_getNodeRanksN = globalenv['singleNode.getNodeRanksN']
G_r = ListVector(G)  # dict to named vector
adj_mat = np.array(adj_mat.data)
ranks = {}

# Walk through all the nodes in S module
logging.debug('Get the single-node encoding node ranks starting from each node.')

for node in S_perturbed_nodes:
    ind = [i+1 for i in range(0, len(G_r.names)) if G_r.names[i] == node][0]
    # Probability diffusion starting from node with index ind
    ranks[node] = singleNode_getNodeRanksN(
        n=ind,
        G=G_r,
        p1=1.0,
        thresholdDiff=0.01,
        adj_mat=adj_mat,
        S=S_perturbed_nodes,
        num_misses=np.log2(len(G_r)),
        verbose=True)

# Vector ranks contains encodings for each node in S_perturbed_nodes
# Prepare for calling R function
mle_getPtBSbyK = globalenv['mle.getPtBSbyK']
ranks_r = ListVector(ranks)  # dict to named vector
S_perturbed_nodes_r = StrVector(S_perturbed_nodes)  # mimics unlist method

# Convert to bitstring
# Get the bitstring associated with the disease module's metabolites
ptBSbyK = mle_getPtBSbyK(S_perturbed_nodes_r, ranks_r)

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
    data_mx_pvals = np.zeros((1,1))

# If we have here specific Patient ID the function will calculate
# Fisher fishers.Info and varPvalue but we'll pass NULL for now
try:
    ptID = data_mx_pvals.columns[0]
except AttributeError:
    ptID = r('NULL')


# Prepare for calling R function
mle_getEncodingLength = globalenv['mle.getEncodingLength']
res = mle_getEncodingLength(ptBSbyK, data_mx_pvals.T, r('NULL'), G_r)
# Returns a subset of nodes that are highly connected
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
ind_F = ind_F.index.values
F_info = res.loc[ind_F, :]

# You can interpret the probability assigned to this metabolite set by
# comparing it to a null encoding algorithm, which uses fixed-length codes
# for all metabolites in the set. The "d.score" is the difference in bitlength
# Significance theorem, we can estimate the upper bounds on a p-value by
# 2^-d.score.
p_value_F = 2.0**(-F_info.iloc[0]['d.score'])

# Prepare for calling R function
get_node_names = globalenv['getNodeNames']

# All metabolites in S
print(S_perturbed_nodes)
# All metabolites in the bitstring
print(get_node_names(ptBSbyK, int(ind_F)))

# Just the F metabolites that are in S_arg that were were "found"

# F_arr = ptBSbyK.rx2(int(ind_F)) <- np.array of encodings
# Fs = np.where(F_arr == 1) <- np.array of indices

# WORKAROUND
Fs = list(get_node_names(ptBSbyK, int(ind_F), found=True))

logging.debug('Set of highly-connected perturbed metabolites F = {} with p-value = {}'.format(Fs, p_value_F))

kmcm_probability = 2.0**(-len(F_info.iloc[0]['optimalBS']))
optimal_bitstring = F_info.iloc[0]['optimalBS']

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
