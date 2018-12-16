#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the adaptive network walk
#' made by the diffusion algorithm trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param data_mx - The matrix that gives the perturbation strength (z-score) for all variables (columns) for each patient (rows).
#' @param ptID - The rowname in pvals associated with the patient being processed.
#' @param perms - The list of permutations calculated over all possible starting nodes, across all metabolites in data.
#' @param kmx - The maximum size of variable sets for which you want to calculate probabilities.
#' @export mle.getPtBSbyK
#' @examples
#' # Read in any network via its adjacency matrix
#' tmp = as.matrix(read.table("adjacency_matrix.txt", sep="\t", header=TRUE))
#' colnames(tmp) = rownames(tmp)
#' ig = graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))  # Must have this declared as a GLOBAL variable!!!!!
#' # Set other tuning parameters
#' p0=0.1  # 10% of probability distributed uniformly
#' p1=0.9  # 90% of probability diffused based on edge weights in networks
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = names(V(ig)$name)
#' # Get node permutations for graph
#' perms = list()
#' for (n in 1:length(G)) {
#'     print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
#'     perms[[names(G)[n]]] = mle.getPermN(n, G)
#' }
#' # Decide what the largest subset size you will consider will be
#' kmx = 20
#' # Load your patient data (p features as rows x n observations as columns)
#' # data_mx = read.table("/your/own/data.txt", sep="\t", header=TRUE)
#' data(testData)
#' data_mx = testData
#' # Get bitstrings associated with each patient's top kmx variable subsets
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'     ptID = colnames(data_mx)[pt]
#'     ptBSbyK[[ptID]] = mle.getPtBSbyK(data_mx, ptID, perms, kmx)
#' }
mle.getPtBSbyK = function(data_mx, ptID, perms, kmx) {
  pt.BSbyK = list()
  pt.sig.nodes = rownames(data_mx)[order(abs(data_mx[,ptID]), decreasing = TRUE)][1:kmx]
  pt.byK = list()
  for (k in 2:kmx) {
    sig.nodes = pt.sig.nodes[1:k]
    pt.bitString = list()
    for (p in 1:length(sig.nodes)) {
      pt.bitString[[sig.nodes[p]]] = as.numeric(perms[[sig.nodes[p]]] %in% sig.nodes)
      names(pt.bitString[[sig.nodes[p]]]) = perms[[sig.nodes[p]]]
      ind = which(pt.bitString[[sig.nodes[p]]]==1)
      pt.bitString[[sig.nodes[p]]] = pt.bitString[[sig.nodes[p]]][1:ind[length(ind)]]
    }
    bestInd = vector("numeric", length(sig.nodes))
    for (p in 1:length(sig.nodes)) {
      bestInd[p] = sum(which(pt.bitString[[p]]==1))
    }
    pt.byK[[k-1]] = pt.bitString[[which.min(bestInd)]]
  }
  pt.BSbyK[[ptID]] = pt.byK
  return(pt.BSbyK)
}

