#!/usr/bin/env Rscript --vanilla
library(methods)
library(argparser)
require(huge)
require(CTD)
require(MASS)

p <- arg_parser("Connect The Dots - Find most connected sub-graph from the set of graphs")
# Add a positional argument
p <- add_argument(p, "--experimental", help="Experimental dataset file name", default = 'vignettes/arg.csv')
p <- add_argument(p, "--control", help="Control dataset file name", default = 'vignettes/control.csv')
# Add a flag
p <- add_argument(p, "--column_name", help="Name of the column containing concentrations")
p <- add_argument(p, "--kmx", help="Number of highly perturbed nodes to consider", default=15)

p <- add_argument(p, "--out", help="output file name")
argv <- parse_args(p)


## I.II: Learn a graph from data.
#Note all code chunks in sections I - IV may rely on lines in previous code
#chunks, so do not empty your environment between code chunks.
#```{r learn_graph}
###data(Miller2015)  # TODO: Add support for user's input file instead od this
# Only include metabolites that are present in >90% reference samples.
###fil.rate=as.numeric(Miller2015$`Times identifed in all 200 samples`[-1])/200
###names(fil.rate) = rownames(Miller2015)[-1]
###data_mx = Miller2015[,grep("IEM_", colnames(Miller2015))]
###data_mx = data_mx[which(fil.rate>0.90), ]
###dim(data_mx)
# Remove any metabolites where any profile has a z-score > 1000.
# These are likely imputed raw values that were not z-scored.
###rmMets = names(which(apply(data_mx, 1, function(i) any(i>20))))
###if (length(rmMets)>0) {
###  data_mx = data_mx[-which(rownames(data_mx) %in% rmMets),]
###}
###dim(data_mx)

# Get data from all patients with Argininemia
###diags = Miller2015["diagnosis", grep("IEM", colnames(Miller2015))]
###experimental_df = data_mx[,which(diags=="Argininemia")]
# TODO: Align these two 196 vs 346 observations!!!!
experimental_df <- read.csv(file = argv$experimental, check.names=FALSE, row.names=1)
experimental_df[] <- lapply(experimental_df, as.character)  # TODO remove?
# Add surrogate disease and surrogate reference profiles based on 1 standard
# deviation around profiles from real patients to improve rank of matrix when
# learning Gaussian Markov Random Field network on data. While surrogate
# profiles are not required, they tend to learn less complex networks
# (i.e., networks with less edges) and in faster time.
control_data <- read.csv(file = argv$control, check.names=FALSE, row.names=1)

target_patients = colnames(experimental_df)
ind = which(diags=="No biochemical genetic diagnosis")
experimental_df=data.surrogateProfiles(experimental_df, 1, ref_data = control_data) #data_mx[,ind])
dim(experimental_df)



# Write graph to file
###write.matrix(arg_icov, file = 'arg_icov.csv', sep = '\t')
arg_icov2 <- read.csv(file = 'arg_icov.csv', sep = '\t', check.names=FALSE) # TODO: Check what check names
rownames(arg_icov2) = colnames(arg_icov2)
arg_icov2 = as.matrix(arg_icov2)
# Convert adjacency matrices to an igraph object.
ig_arg = graph.adjacency(arg_icov2, mode="undirected", weighted=TRUE,
                         add.colnames="name")
#print(ig_arg)


# IV. The Encoding Process
#We're going to go back to our data using the Arginase deficiency network
#model, and the Miller et al (2015) dataset.
## IV.0 Re-define the Arginase deficiency network model

print(ig_arg)
adj_mat = as.matrix(get.adjacency(ig_arg, attr="weight"))
G=vector(mode="list", length=length(V(ig_arg)$name))
G[1:length(G)] = 0
names(G) = V(ig_arg)$name


## IV.I Choose your node subset.


# Maximum subset size to inspect
kmx=argv$kmx
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

S_disease_module_ind = which(occurances$Freq >= (length(target_patients) / 2.))
S_disease_module = as.list(as.character(occurances[S_disease_module_ind, "vec_s"]))

# TODO: Replace names(S_arg) with S_disease_module

## IV.II Get k node permutations.
# Get the single-node encoding node ranks starting from each node in the subset
# S_arg.
ranks = list()
for (n in 1:length(S_arg)) {
  ind = which(names(G)==names(S_arg)[n])
  # probability_difussion starting from node with index ind
  ranks[[n]]=singleNode.getNodeRanksN(ind,G,p1=1.0,thresholdDiff=0.01,
                                      adj_mat,S=names(S_arg),
                                      num.misses=log2(length(G)),TRUE)
}
names(ranks) = names(S_arg)
# Vector ranks contains encodings for each node in S_arg
# TODO: Why ranks contain 87 rows instead of 15?!

## IV.III Convert to bitstrings.

# Get the bitstrings associated with the patient's perturbed metabolites in
# "S_arg" based on the node ranks calculated in the previous step, "ranks".
ptBSbyK = mle.getPtBSbyK(names(S_arg), ranks)


## IV.IV Get encoding length of minimum length codeword.
# experimental_df is dataframe with diseases (and surrogates)
# and z-values for each metabolite
ind = which(colnames(experimental_df) %in% names(diags))
data_mx.pvals=apply(experimental_df[,ind], c(1,2),
                    function(i) 2*pnorm(abs(i), lower.tail=FALSE))
ptID = "IEM_1006"
res = mle.getEncodingLength(ptBSbyK, t(data_mx.pvals), ptID, G)
# returns a subset of nosed that are highly connected
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
names(S_arg)
# Which metabolites were in the 8 metabolite subset of patient IEM_1006's
# top 15 perturbed metabolites that had the above p-value?
ptBSbyK[[ind.mx]] # all metabolites in the bitstring
# just the F metabolites that are in S_arg that were were "found"
F = names(which(ptBSbyK[[ind.mx]]==1))
print(F)
print(p_value_F)
