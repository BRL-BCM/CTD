#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the adaptive network walk
#' made by the diffusion algorithm trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param sig.nodes - A character vector of node names.
#' @param perms - The list of permutations calculated over all possible nodes, starting with each node in subset of interest.
#' @return pt.byK - a list of bitstrings, with the names of the list elements the node names of the encoded nodes
#' @export kraft.getPtBSbyK
kraft.getPtBSbyK = function(pt.sig.nodes, perms) {
  pt.byK = list()
  for (k in 2:kmx) {
    sig.nodes = pt.sig.nodes[1:k]
    pt.bitString = list()
    for (p in 1:length(sig.nodes)) {
      pt.bitString[[sig.nodes[p]]] = as.numeric(perms[[sig.nodes[p]]] %in% sig.nodes)
      names(pt.bitString[[sig.nodes[p]]]) = perms[[sig.nodes[p]]]
      ind = which(pt.bitString[[sig.nodes[p]]] == 1)
      pt.bitString[[sig.nodes[p]]] = pt.bitString[[sig.nodes[p]]][1:ind[length(ind)]]
    }
    # Which found the most nodes
    bestInd = vector("numeric", length(sig.nodes))
    for (p in 1:length(sig.nodes)) {
      bestInd[p] = sum(pt.bitString[[p]])
    }
    pt.bitString = pt.bitString[which(bestInd==max(bestInd))]
    # Which found the most nodes soonest
    bestInd = vector("numeric", length(pt.bitString))
    for (p in 1:length(pt.bitString)) {
      bestInd[p] = sum(which(pt.bitString[[p]] == 1))
    }
    pt.byK[[k - 1]] = pt.bitString[[which.max(bestInd)]]
  }
  return(pt.byK)
}

