#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the adaptive network walk
#' made by the diffusion algorithm trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param data - The matrix that gives the perturbation strength (z-score) for all variables (columns) for each patient (rows).
#' @param ptID - The rowname in pvals associated with the patient being processed.
#' @param perms - The list of permutations calculated over all possible starting nodes, across all metabolites in data.
#' @param kmx - The maximum size of variable sets for which you want to calculate probabilities.
#' @export mle.getPtBSbyK
#' @examples
#' mle.getPtBSbyK(data, ptID, permutationByStartNode, kmax)
mle.getPtBSbyK = function(data, ptID, perms, kmx) {
  pt.BSbyK = list()
  pt.sig.nodes = rownames(data)[order(abs(data[,ptID]), decreasing = TRUE)][1:kmx]
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

