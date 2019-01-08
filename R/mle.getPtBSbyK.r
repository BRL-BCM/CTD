#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the adaptive network walk
#' made by the diffusion algorithm trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param data_mx - The matrix that gives the perturbation strength (z-score) for all variables (rows) for each patient (columns).
#' @param ptID - The identifier associated with the patient being processed.
#' @param perms - The list of permutations calculated over all possible nodes, starting with each node in subset of interest.
#' @param kmx - The maximum size of variable sets for which you want to calculate probabilities.
#' @return pt.byK - a list of bitstrings, with the names of the list elements the node names of the encoded nodes
#' @export mle.getPtBSbyK
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Get bitstrings associated with each patient's top kmx variable subsets
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'   ptID = colnames(data_mx)[pt]
#'   ptBSbyK[[ptID]] = mle.getPtBSbyK(data_mx, ptID, perms, kmx)
#' }
mle.getPtBSbyK = function(data_mx, ptID, perms, kmx) {
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
  return(pt.byK)
}

