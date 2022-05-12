#!/usr/bin/env Rscript --vanilla
library(argparser)
require(huge)
require(MASS)
library(rjson)
library(stringr)
library(fs)
require(igraph)
source("./R/mle.getEncodingLength.r")
source("./R/mle.getPtBSbyK.r")
source("./R/data.surrogateProfiles.r")
source("./R/data.imputeData.r")
source("./R/singleNode.getNodeRanksN.r")
source("./R/graph.diffuseP1.r")
source("./R/graph.connectToExt.r")
source("./R/stat.fishersMethod.r")

p <- arg_parser("Connect The Dots - Find the most connected sub-graph")
p <- add_argument(p, "--experimental", help="Experimental dataset file name", default = '')  # data/example_argininemia/experimental.csv
p <- add_argument(p, "--control", help="Control dataset file name", default = '')            # data/example_argininemia/control.csv
p <- add_argument(p, "--adj_matrix", help="CSV with adjacency matrix", default = '')         # data/example_argininemia/adj.csv
p <- add_argument(p, "--disease_module", help="Comma-separated list of graph G nodes to consider when searching for the most connected sub-graph")
p <- add_argument(p, "--kmx", help="Number of highly perturbed nodes to consider. Ignored if disease_module is given.", default=15)
p <- add_argument(p, "--present_in_perc", help="Percentage of patients having metabolite. Ignored if disease_module is given.", default=0.5)
p <- add_argument(p, "--output_name", help="Name of the output JSON file.", default=NA)
p <- add_argument(p, "--out_graph_name", help="Name of the output graph adjecancy CSV file.", default=NA)

argv <- parse_args(p)

# Read input dataframe with experimental (positive, disease) samples
if (file.exists(argv$experimental)){
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
  experimental_df = as.data.frame(experimental_df)
} else{
  target_patients = NA
}


# Read input graph (adjacency matrix)
if (file.exists(argv$adj_matrix)){
  adj_df <- read.csv(file = argv$adj_matrix, check.names=FALSE)
  rownames(adj_df) = colnames(adj_df)
  adj_df = as.matrix(adj_df)
} else{
  experimental = huge(t(experimental_df), method="glasso")
  # This will take several minutes. For a faster option, you can use the
  # "ebic" criterion instead of "stars", but we recommend "stars".
  experimental.select = huge.select(experimental, criterion="stars")
  adj_df = as.matrix(experimental.select$opt.icov)
  diag(adj_df) = 0
  rownames(adj_df) = rownames(experimental_df)
  colnames(adj_df) = rownames(experimental_df)
  if (!is.na(argv$out_graph_name)){
    write.table(adj_df, file = argv$out_graph_name, row.names=FALSE, sep=',')
  }
}
# Convert adjacency matrices to an igraph object.
igraph = graph.adjacency(adj_df, mode="undirected", weighted=TRUE,
                         add.colnames="name")
# The Encoding Process
adj_mat = as.matrix(get.adjacency(igraph, attr="weight"))
G=vector(mode="list", length=length(V(igraph)$name))
G[1:length(G)] = 0
names(G) = V(igraph)$name

## Choose node subset
kmx=argv$kmx  # Maximum subset size to inspect
S_set = list()

if (is.na(argv$disease_module)){
  for (pt in target_patients) {
    sel = experimental_df[[pt]]
    temp = experimental_df[order(sel, decreasing=TRUE), ]
    S_patient=rownames(temp)[1:kmx]
    S_set<-append(S_set, S_patient)
  } # Created list containing top kmx metabolites for every taget user
  vec_s = unlist(S_set)
  occurances = as.data.frame(table(vec_s))
  # Keep in the disease_module the metabolites perturbed in at least 50% patients
  S_disease_module_ind = which(occurances$Freq >= (length(target_patients) * argv$present_in_perc))
  S_disease_module = as.list(as.character(occurances[S_disease_module_ind, "vec_s"]))
} else
{
  S_disease_module = as.list(unlist(str_split(argv$disease_module, ',')))
  # TODO: Check if all given nodes are in G!
}
## Get k node permutations
# Get the single-node encoding node ranks starting from each node in the subset
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

## Convert to bitstrings
# Get the bitstrings associated with the disease module's metabolites
ptBSbyK = mle.getPtBSbyK(unlist(S_disease_module), ranks)

## Get encoding length of minimum length code word.
# experimental_df is dataframe with diseases (and surrogates)
# and z-values for each metabolite
# TODO: If graph and disease module are given do we still need experimental_df?
if (file.exists(argv$experimental)){
  ind = which(colnames(experimental_df) %in% target_patients)
  data_mx.pvals=apply(experimental_df[,ind], c(1,2),
                    function(i) 2*pnorm(abs(i), lower.tail=FALSE))
  # p-value is area under curve of normal distribution on the right of the
  # specified z-score. pnorm generates normal distribution with mean=0, std=1,
  # which is exactly what z-scores are
}else{
  data_mx.pvals = list(0)
}
ptID = colnames(data_mx.pvals)[1] # If we have here specific Patient ID the function will calculate
            # Fisher fishers.Info and varPvalue
res = mle.getEncodingLength(ptBSbyK, t(data_mx.pvals), ptID, G)
# returns a subset of nodes that are highly connected
ind.mx = which.max(res$d.score)
ind.mx = which(res$d.score == max(res$d.score) )
highest_dscore_paths = res[ind.mx,]

## Locate encoding (F) with best d-score
# Tiebraker 1: If several results have the same d-score take one with longest BS
a1 = nchar(highest_dscore_paths$optimalBS)
a2 = max(nchar(highest_dscore_paths$optimalBS))
ind_F = which(nchar(highest_dscore_paths$optimalBS) == 
                max(nchar(highest_dscore_paths$optimalBS)))
ind_F = highest_dscore_paths[ind_F,]
# Tiebraker 2: Take the one with the largest subsetSize
index_highest = which(ind_F$subsetSize == 
                        max(ind_F$subsetSize))
ind_F = ind_F[index_highest,]
ind_F = rownames(ind_F)
F_info = res[ind_F,]

# You can interpret the probability assigned to this metabolite set by
# comparing it to a null encoding algorithm, which uses fixed-length codes
# for all metabolites in the set. The "d.score" is the difference in bitlength
# Significance theorem, we can estimate the upper bounds on a p-value by
# 2^-d.score.

p_value_F = 2^-F_info$d.score
S_disease_module # All metabolites in S
ptBSbyK[[ind_F]] # all metabolites in the bitstring
# just the F metabolites that are in S_arg that were were "found"
F_arr = ptBSbyK[[strtoi(ind_F)]]
F = which(F_arr==1)
F = names(F)
print(paste('Set of highly-connected perturbed metabolites F = {', toString(F), 
            '} with p-value = ', p_value_F))

out_dict <- list(most_connected_nodes = F,p_value = p_value_F)
res_json = toJSON(out_dict, indent = 4)
if (is.na(argv$output_name)){
  if (argv$experimental == ''){
    outfname = fs::path_file(argv$adj_matrix)
  }else{
    outfname = fs::path_file(argv$experimental)
  }
  outfname = str_replace(outfname, 'csv', 'json')
}else{
  outfname = argv$output_name
}
write(res_json, outfname)