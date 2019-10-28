#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the adaptive network walk
#' made by the diffusion algorithm trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param S - A character vector of node names describing the node subset to be encoded.
#' @param perms - The list of permutations calculated over all possible nodes, starting with each node in subset of interest.
#' @return pt.byK - a list of bitstrings, with the names of the list elements the node names of the encoded nodes
#' @export mle.getPtBSbyK_memory
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Get bitstrings associated with each patient's top kmx variable subsets
#' kmx = 15
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'   S = data_mx[order(abs(data_mx[,pt]), decreasing=TRUE),pt][1:kmx]
#'   ptBSbyK[[ptID]] = mle.getPtBSbyK_memory(S, perms)
#' }
#mle.getPtBSbyK_memory = function(S, perms) {
#  pt.byK = list()
#  for (k in 1:length(S)) {
#    sig.nodes = S[1:k]
#    sig.nodes = sapply(sig.nodes, trimws)
#    pt.bitString = list()
#    for (p in 1:length(sig.nodes)) {
#      pt.bitString[[sig.nodes[p]]] = as.numeric(perms[[sig.nodes[p]]] %in% sig.nodes)
#      names(pt.bitString[[sig.nodes[p]]]) = perms[[sig.nodes[p]]]
#      ind = which(pt.bitString[[sig.nodes[p]]] == 1)
#      pt.bitString[[sig.nodes[p]]] = pt.bitString[[sig.nodes[p]]][1:ind[length(ind)]]
#    }
    # Which found the most nodes
#    bestInd = vector("numeric", length(sig.nodes))
#    for (p in 1:length(sig.nodes)) {
#      bestInd[p] = sum(pt.bitString[[p]])
#    }
#    pt.bitString = pt.bitString[which(bestInd==max(bestInd))]
    # Which found the most nodes soonest
#    bestInd = vector("numeric", length(pt.bitString))
#    for (p in 1:length(pt.bitString)) {
#      bestInd[p] = sum(which(pt.bitString[[p]] == 1))
#          }
#    pt.byK[[k]] = pt.bitString[[which.min(bestInd)]]
#  }
#  return(pt.byK)
#}


mle.getPtBSbyK_memory = function(S, perms) {
  perms2 = perms[which(names(perms) %in% S)]
  pt.bitString = list()
  for (p in 1:length(perms2)) {
    pt.bitString[[S[p]]] = as.numeric(perms2[[S[p]]] %in% S)
    names(pt.bitString[[S[p]]]) = perms2[[S[p]]]
    ind = which(pt.bitString[[S[p]]] == 1)
    pt.bitString[[S[p]]] = pt.bitString[[S[p]]][1:ind[length(ind)]]
  }
  # For each k, find the perms that found at least k in S. Which found the first k soonest?
  pt.byK = list()
  for (k in 1:length(S)) {
    pt.byK_tmp = pt.bitString
    # Which found at least k nodes
    bestInd = vector("numeric", length(S))
    for (p in 1:length(S)) {
      bestInd[p] = sum(pt.bitString[[p]])
    }
    if (length(which(bestInd>=k))>0) {
      pt.byK_tmp = pt.byK_tmp[which(bestInd>=k)]
      # Which found the most nodes soonest
      bestInd = vector("numeric", length(pt.byK_tmp))
      for (p in 1:length(pt.byK_tmp)) {
        bestInd[p] = sum(which(pt.byK_tmp[[p]] == 1)[1:k])
      }
      maxInd = max(which(pt.byK_tmp[[which.min(bestInd)]] == 1)[1:k])
      pt.byK[[k]] = pt.byK_tmp[[which.min(bestInd)]][1:maxInd]
    } else {
      pt.byK[[k]] = pt.byK_tmp[[which.max(bestInd)]]
    }
  }

  return(pt.byK)
}
