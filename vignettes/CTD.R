#!/usr/bin/env Rscript --vanilla
library(methods)
library(argparser)
require(huge)
require(CTD)
require(MASS)

p <- arg_parser("Connect The Dots - Find the most connected sub-graph from the set of graphs")
# Add a positional argument
p <- add_argument(p, "--experimental", help="Experimental dataset file name", default = 'data/example_argininemia/arg.csv')
p <- add_argument(p, "--control", help="Control dataset file name", default = 'data/example_argininemia/control.csv')
p <- add_argument(p, "--adj_matrix", help="CSV with adjecancy matric", default = 'data/example_argininemia/arg_adj.csv')
# Add a flag
p <- add_argument(p, "--column_name", help="Name of the column containing concentrations")
p <- add_argument(p, "--kmx", help="Number of highly perturbed nodes to consider", default=15)
p <- add_argument(p, "--out", help="output file name")
argv <- parse_args(p)

# Read input dataframe with experimental (positive, disease) samples
experimental_df <- read.csv(file = argv$experimental, check.names=FALSE, row.names=1)
experimental_df[] <- lapply(experimental_df, as.character)  # TODO remove?

# Read input dataframe with control samples
control_data <- read.csv(file = argv$control, check.names=FALSE, row.names=1)

target_patients = colnames(experimental_df)

# TODO: Do we need surrogates?
# Add surrogate disease and surrogate reference profiles based on 1 standard
# deviation around profiles from real patients to improve rank of matrix when
# learning Gaussian Markov Random Field network on data. While surrogate
# profiles are not required, they tend to learn less complex networks
# (i.e., networks with less edges) and in faster time.
experimental_df=data.surrogateProfiles(experimental_df, 1, ref_data = control_data)
dim(experimental_df)

# Read input graph (adjacency matrix)
adj_df <- read.csv(file = argv$adj_matrix, sep = '\t', check.names=FALSE)
rownames(adj_df) = colnames(adj_df)
adj_df = as.matrix(adj_df)
# Convert adjacency matrices to an igraph object.
igraph = graph.adjacency(adj_df, mode="undirected", weighted=TRUE,
                         add.colnames="name")

# IV. The Encoding Process
adj_mat = as.matrix(get.adjacency(igraph, attr="weight"))
G=vector(mode="list", length=length(V(igraph)$name))
G[1:length(G)] = 0
names(G) = V(igraph)$name

## IV.I Choose your node subset.
kmx=argv$kmx  # Maximum subset size to inspect
# Get our node subset associated with the $KMX highest perturbed (up or down)
# in our first Arginase deficiency sample.
S_set = list()

experimental_df = as.data.frame(experimental_df)
for (pt in target_patients) {
  sel = experimental_df[[pt]]
  temp = experimental_df[order(sel), ]
  S_patient=rownames(temp)[1:kmx]
  S_set<-append(S_set, S_patient)
}
# TODO: Take union >50%
vec_s = unlist(S_set)
occurances = as.data.frame(table(vec_s))

S_disease_module_ind = which(occurances$Freq >= (length(target_patients) / 2.))  # TODO: Create apearance_threshold variable
S_disease_module = as.list(as.character(occurances[S_disease_module_ind, "vec_s"]))

## IV.II Get k node permutations.
# Get the single-node encoding node ranks starting from each node in the subset
# S_arg.
ranks = list()
for (n in 1:length(S_disease_module)) {
  ind = which(names(G)==S_disease_module[n])
  # probability_difussion starting from node with index ind
  ranks[[n]]=singleNode.getNodeRanksN(ind,G,p1=1.0,thresholdDiff=0.01,
                                      adj_mat,S=S_disease_module,
                                      num.misses=log2(length(G)),TRUE)
}
names(ranks) = S_disease_module
# Vector ranks contains encodings for each node in S_disease_module

## IV.III Convert to bitstrings.
# Get the bitstrings associated with the disease module's metabolites
ptBSbyK = mle.getPtBSbyK(unlist(S_disease_module), ranks)

## IV.IV Get encoding length of minimum length code word.
# experimental_df is dataframe with diseases (and surrogates)
# and z-values for each metabolite
ind = which(colnames(experimental_df) %in% names(diags))
data_mx.pvals=apply(experimental_df[,ind], c(1,2),
                    function(i) 2*pnorm(abs(i), lower.tail=FALSE))

ptID = "IEM_1006" # TODO: Why we still need this?!
res = mle.getEncodingLength(ptBSbyK, t(data_mx.pvals), ptID, G)
# returns a subset of nodes that are highly connected
ind.mx = which.max(res$d.score)
res[ind.mx,]

## IV.V Get p-value of variable length encoding vs. fixed length encoding.

# You can interpret the probability assigned to this metabolite set by
# comparing it to a null encoding algorithm, which uses fixed-length codes
# for all metabolites in the set. The "d.score" is the difference in bitlength
# between the null and alternative encoding models. Using the Algorithmic
# Significance theorem, we can estimate the upper bounds on a p-value by
# 2^-d.score.
p_value_F = 2^-res[ind.mx,"d.score"]
# All metabolites in S_arg
S_disease_module
# Which metabolites were in the 8 metabolite subset of patient IEM_1006's
# top 15 perturbed metabolites that had the above p-value?
ptBSbyK[[ind.mx]] # all metabolites in the bitstring
# just the F metabolites that are in S_arg that were were "found"
F = names(which(ptBSbyK[[ind.mx]]==1))
print(F)
print(p_value_F)
