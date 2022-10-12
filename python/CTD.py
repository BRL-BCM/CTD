import argparse
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

huge = importr("huge")

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

print(argv.adj_matrix)
# Read input dataframe with experimental (positive, disease) samples
# if (file.exists(argv.experimental)){
#   experimental_df <- read.csv(file = argv$experimental, check.names=FALSE, row.names=1)
#   experimental_df[] <- lapply(experimental_df, as.character)  # TODO remove?
#   # Read input dataframe with control samples
#   control_data <- read.csv(file = argv$control, check.names=FALSE, row.names=1)
#   target_patients = colnames(experimental_df)