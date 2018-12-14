require(MolPhenoMatch)

# 1. DATA LOAD ---------------------------------------------------------------
# 1.1 Load the raw data into R (N by p matrix, where N are the samples of p features). Please be advised
# that our graphical learning phase works best with data or fitted data that are normally distributed.
# We use Gaussian Markov Random Field models to approximate the partial correlation structure of the p features.
# :: Data matrix::
#     Rows: patients/samples
#     Columns: features/variables.
# data.raw = read.table("/path/to/data/matrix.txt", sep="\t", header=TRUE)
data(testData)

# 1.2 If raw data are already fit to a normal distribution, skip step 1.2 below.
isControlSample = c(rep(0,nrow(testData)/2), rep(1,nrow(testData)/2))
data.zscore = data.fitDataToNormal(testData,isControlSample)

# 1.3 Convert the z-scores to p-values, tail agnostic. Make sure rows are patients/samples, and columns are features/variables.
data.pvals = apply(data.zscore, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
data.pvals = t(data.pvals)
dim(data.pvals)






# 2. EXPERIMENT PREP ---------------------------------------------------------
# 2.1 Learn your background knowledge graph on your data, using a Gaussian Markov Random Field Models.
# 2.1.1 For N > p, learn the inverse of the covariante matrix directly from the data.
require(MASS)
covmat = cov(t(data.zscore))
inv_covmat = ginv(covmat)
diag(inv_covmat) = 0;
colnames(inv_covmat) = rownames(data.zscore)
ig = graph.adjacency(as.matrix(inv_covmat), mode="undirected", weighted=TRUE, add.colnames='name')


# 2.1.2 For N << p data sources, we recommend using "surrogate" profiles to improve the rank of N relative to p.
# Then, learn the inverse of the covariance matrix directly.



require(MASS)
covmat = cov(t(data.zscore))
inv_covmat = ginv(covmat)
diag(inv_covmat) = 0;
colnames(inv_covmat) = rownames(data.zscore)
ig = graph.adjacency(as.matrix(inv_covmat), mode="undirected", weighted=TRUE, add.colnames='name')




# 2.1.3 For N < p, the inverse of the covariance matrix can be approximated using the graphical LASSO.
# We recommend setting the regularization parameter, lambda, to a very small number, to preserve the majority
# of information embedded in the partial correlation structure.
# Features should be columns, patients should be rows.
require(huge)
inv_covmatt = huge(t(data.zscore), method="glasso", lambda = 0.01)
inv_covmat = as.matrix(inv_covmatt$icov[[1]])
diag(inv_covmat) = 0;
colnames(inv_covmat) = rownames(data.zscore)
ig = graph.adjacency(as.matrix(inv_covmat), mode="undirected", weighted=TRUE, add.colnames='name')


# 2.2 Set your global tuning parameter variables: p0, p1, thresholdDiff, and thresholdDrawT.
p0=0.1  # This is the percentage of probability that gets divided uniformly between the variables in the background knowledge graph ("ig") that was learned in 2.1.
p1=0.9  # This is the percentage of probability that gets divided preferentially (depending on the connectivity) between the variables in the background knowledge graph that was learned in 2.1.
thresholdDiff=0.01  # This is the threshold of probability for which the probability diffusion algorithm stops recursively diffusing.
kmax=5  # This is the size of the largest subset you will encode against your background knowledge graph.
igraphObjectG = vector(mode="list", length=length(V(ig)$name))  # Global placeholder for the probabilities that are assigned to each variable in the background graph.
names(igraphObjectG) = V(ig)$name
adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))  # Global copy of the background knowledge graph's adjacency matrix.

# 2.3 It is strongly advised that you save your progress up until this point into a smartly named RData file, so you can begin at stage 3 at any given time.
path= paste(getwd(), "/projectFolder", sep="")
save.image(sprintf("%s/experimentPrep.RData", path))




# PRE-COMPUTE NODE PERMUTATIONS UP FRONT --------------------------------------------------------
# If you have use of a computing cluster, each permutation can be run in parallel and then collated
# into one .RData file that contains the list object "permutationByStartNode" If you do not have a computing cluster,
# and in particular, if the number of nodes in your graph is manageable, running in serial as shown
# below is a safe option.
permutationByStartNode = list()
for (n in 1:length(igraphObjectG)) {
  all_nodes = names(igraphObjectG)
  current_node_set = NULL
  stopIterating=FALSE
  startNode = all_nodes[n]
  currentGraph = igraphObjectG
  while (stopIterating==FALSE) {
    current_node_set = c(current_node_set, startNode)
    #STEP 4: Diffuse p0 and p1 to connected nodes from current draw, startNode node
    #For unseen nodes, clear probabilities and add probability (p0/#unseen nodes)
    for (t in 1:length(currentGraph)) {
      currentGraph[[t]] = 0 #set probabilities of all nodes to 0
    }
    #determine base p0 probability
    baseP = p0/(length(currentGraph)-length(current_node_set))
    for (t in 1:length(currentGraph)) {
      if (!(names(currentGraph[t]) %in% current_node_set)) {
        currentGraph[[t]] = baseP  #set probabilities of unseen nodes to diffused p0 value, baseP
      } else {
        currentGraph[[t]] = 0
      }
    }
    p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
    currentGraph = graph.diffuseP1(p1, startNode, currentGraph, currentGraph[current_node_set], 1, verbose=FALSE)
    p1_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
    if (abs(p1_event-1)>thresholdDiff) {
      extra.prob.to.diffuse = 1-p1_event
      currentGraph[names(current_node_set)] = 0
      currentGraph[!(names(currentGraph) %in% names(current_node_set))] = unlist(currentGraph[!(names(currentGraph) %in% names(current_node_set))]) + extra.prob.to.diffuse/sum(!(names(currentGraph) %in% names(current_node_set)))
    }
    #Set startNode to a node that is the max probability in the new currentGraph
    maxProb = names(which.max(currentGraph))
    # Break ties: When there are ties, choose amongst winners for highest degree.
    # Just pick top of alphabet in maxProb vector.
    startNode = names(currentGraph[maxProb[1]])
    if (length(c(startNode,current_node_set))>=(length(V(igraphTestG)$name))) {
      current_node_set = c(current_node_set, startNode)
      stopIterating = TRUE
    }
  }
  permutationByStartNode[[all_nodes[n]]] = current_node_set
}



# GET PATIENT BITSTRINGS --------------------------------------------------------
# If data load (Section 1) and experiment prep (Section 2) are already saved, load environment here.
# Specify a path where patient-specific data will be stored.
# setwd("/path/to/projectDir")
setwd(paste(getwd(), "/projectFolder", sep=""))
load("experimentPrep.RData")
# Takes about 2-3 minutes for test dataset
for (patient in 1:nrow(data.pvals)) {
  mle.getPatientBitstrings(ig, path, data.pvals, patient)
}



# GET MINIMUM ENCODING LENGTH ---------------------------------------------
# You can look at the particular encoding lengths of a given patient with this function. The mutual information metric in the next section uses this function, too.
# Subsets that were highly connected will have more 1's in the BitString, and larger d-scores (See the "d.entropy" column in the results object, below).
patientNum=rownames(data.pvals)[10]
results = mle.getEncodingLength(ig, path, data.pvals, patientNum)
results = results[order(results[,"d.score"], decreasing = TRUE),]
sizeSubset = results[1,"subsetSize"]



# GET PATIENT SIMILARITY MATRIX -------------------------------------------
# Note: Each comparison takes an average of 2 seconds. There are 50*50 pairwise comparisons total.
#       Expected execution time (if executed in serial) is currently 75-90 minutes for the test dataset.
patientSimilarity = matrix(0, nrow=nrow(data.pvals), ncol=nrow(data.pvals))
rownames(patientSimilarity) = rownames(data.pvals)
colnames(patientSimilarity) = rownames(data.pvals)
for (p1 in 1:nrow(patientSimilarity)) {
  pt1 = rownames(data.pvals)[p1]
  print(pt1)
  for (p2 in 1:ncol(patientSimilarity)) {
    pt2 = rownames(data.pvals)[p2]
    patientSimilarity[pt1,pt2] = mle.getPatientSimilarity(ig, path, data.pvals, pt1, pt2)
  }
}
# Normalize scale by row so min is 0 and max is 1.
for (p1 in 1:nrow(patientSimilarity)) {
  tmp = as.matrix(patientSimilarity[p1,])
  patientSimilarity[p1,] = (tmp-min(tmp))/(max(tmp)-min(tmp))
}
plot.ptSimilarity(patientSimilarity, path="/Users/lillian.rosa/Downloads/Thesis/projectFolder")
