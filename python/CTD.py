import os
import sys
import argparse
import pandas as pd
import numpy as np
from igraph import Graph, _igraph
from rpy2.robjects import r, pandas2ri, packages, globalenv


def eprint(args):
    sys.stderr.write(str(args) + "\n")

huge = packages.importr("huge")
r.source("./CTD/R/data.surrogateProfiles.r")
r.source("./CTD/R/data.imputeData.r")
pandas2ri.activate()

p = argparse.ArgumentParser(description="Connect The Dots - Find the most connected sub-graph")
p.add_argument("--experimental", help="Experimental dataset file name", default = '')  # data/example_argininemia/experimental.csv
p.add_argument("--control", help="Control dataset file name", default = '')            # data/example_argininemia/control.csv
p.add_argument("--adj_matrix", help="CSV with adjacency matrix", default = '')         # data/example_argininemia/adj.csv
p.add_argument("--s_module", help="Comma-separated list or path to CSV of graph G nodes to consider when searching for the most connected sub-graph")
p.add_argument("--kmx", help="Number of highly perturbed nodes to consider. Ignored if S module is given.", default=15)
p.add_argument("--present_in_perc_for_s", help="Percentage of patients having metabolite for selection of S module. Ignored if S module is given.", default=0.5)
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

    data_surrogateProfiles = globalenv['data.surrogateProfiles']
    experimental_df = data_surrogateProfiles(data=experimental_df, std=1, ref_data=control_data)
else:
    target_patients = np.nan

# Read input graph (adjacency matrix)
if os.path.exists(argv.adj_matrix):
    adj_df = pd.read_csv(argv.adj_matrix)
    adj_df.index = adj_df.columns

else:

    experimental = huge.huge(experimental_df.to_numpy().T, method='glasso')
    # This will take several minutes. For a faster option, you can use the
    # "ebic" criterion instead of "stars", but we recommend "stars".
    eprint('Starting huge.select.')
    experimental_select = huge.huge_select(experimental, criterion=argv.glasso_criterion)
    adj_df = [values for (name, values) in experimental_select.items() if name == 'opt.icov'][0]
    np.fill_diagonal(adj_df, 0)
    adj_df = pd.DataFrame(adj_df, index=experimental_df.index, columns=experimental_df.index)
    if argv.out_graph_name:
        adj_df.to_csv(argv.out_graph_name, index=False)


eprint('Convert adjacency matrices to an igraph object.')
try:
    igraph = Graph.Weighted_Adjacency(adj_df, mode='undirected')
except _igraph.InternalError:
    igraph = Graph.Weighted_Adjacency(adj_df, mode='upper')

print(igraph)
