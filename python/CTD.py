import os
import sys
import argparse
import pandas as pd
import numpy as np
from igraph import Graph
from rpy2.robjects import r, pandas2ri, numpy2ri, packages, globalenv, ListVector
from collections import Counter


def eprint(args):
    sys.stderr.write(str(args) + "\n")

huge = packages.importr("huge")
r.source("./CTD/R/data.surrogateProfiles.r")
r.source("./CTD/R/data.imputeData.r")
r.source("./CTD/R/singleNode.getNodeRanksN.r")
r.source("./CTD/R/graph.diffuseP1.r")
r.source("./CTD/R/graph.connectToExt.r")
r.source("./CTD/R/stat.fishersMethod.r")

pandas2ri.activate()

p = argparse.ArgumentParser(description="Connect The Dots - Find the most connected sub-graph")
p.add_argument("--experimental", help="Experimental dataset file name", default = '')  # data/example_argininemia/experimental.csv
p.add_argument("--control", help="Control dataset file name", default = '')            # data/example_argininemia/control.csv
p.add_argument("--adj_matrix", help="CSV with adjacency matrix", default = '')         # data/example_argininemia/adj.csv
p.add_argument("--s_module", help="Comma-separated list or path to CSV of graph G nodes to consider when searching for the most connected sub-graph")
p.add_argument("--kmx", help="Number of highly perturbed nodes to consider. Ignored if S module is given.", default=15, type=int)
p.add_argument("--present_in_perc_for_s", help="Percentage of patients having metabolite for selection of S module. Ignored if S module is given.", default=0.5, type=float)
p.add_argument("--output_name", help="Name of the output JSON file.")
p.add_argument("--out_graph_name", help="Name of the output graph adjecancy CSV file.")
p.add_argument("--glasso_criterion", help="Graph-ical Lasso prediction of the graph criterion. stars is default, ebic is faster.", default='stars')
argv = p.parse_args()

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
    eprint('Starting huge.select.')
    experimental_select = huge.huge_select(experimental, criterion=argv.glasso_criterion)
    adj_df = [values for (name, values) in experimental_select.items() if name == 'opt.icov'][0]
    np.fill_diagonal(adj_df, 0)
    adj_df = pd.DataFrame(adj_df, index=experimental_df.index, columns=experimental_df.index)  # keep as DataFrame for now
    if argv.out_graph_name:
        adj_df.to_csv(argv.out_graph_name, index=False)

eprint('Convert adjacency matrices to an igraph object.')
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

eprint('Selected perturbed nodes, S = {}'.format(S_perturbed_nodes))

# Check if all nodes from the s_module are in graph
for node in S_perturbed_nodes:
    if node not in G:
        eprint('Node "{}" not in graph. Exiting program.'.format(node))
        exit(1)

# Walk through all the nodes in S module
eprint('Get the single-node encoding node ranks starting from each node.')

# Prepare for calling R function
singleNode_getNodeRanksN = globalenv['singleNode.getNodeRanksN']
G_r = ListVector(G)
adj_mat = np.array(adj_mat.data)
ranks = {}

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
ranks_r = ListVector(ranks)

# Convert to bitstring
# Get the bitstring associated with the disease module's metabolites
