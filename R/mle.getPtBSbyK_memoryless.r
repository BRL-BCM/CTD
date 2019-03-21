#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the memoryless network walker
#' trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param S - A character vector of node names describing the node subset to be encoded.
#' @param perms - The list of permutations calculated over all possible nodes, starting with each node in subset of interest.
#' @param num.misses - The number of misses tolerated by the network walker before path truncation occurs.
#' @return pt.byK - a list of bitstrings, with the names of the list elements the node names of the encoded nodes
#' @export mle.getPtBSbyK_memoryless
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Get bitstrings associated with each patient's top kmx variable subsets
#' kmx = 15
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'   S = data_mx[order(abs(data_mx[,pt]), decreasing=TRUE),pt][1:kmx]
#'   ptBSbyK[[ptID]] = mle.getPtBSbyK_memoryless(S, perms)
#' }
mle.getPtBSbyK_memoryless = function(S, perms, num.misses=NULL) {
  if (is.null(num.misses)) {
    num.misses = ceiling(log2(length(perms)))
  }
  pt.byK = list()
  for (k in 1:length(S)) {
    sig.nodes = S[1:k]
    sig.nodes = sapply(sig.nodes, trimws)
    pt.bitString = list()
    for (p in 1:length(sig.nodes)) {
      miss = 0
      for (ii in 1:length(perms[[sig.nodes[p]]])) {
        ind_t = as.numeric(perms[[sig.nodes[p]]][ii] %in% sig.nodes)
        if (ind_t==0) {
          miss = miss + 1
          if (miss >= num.misses) {
            thresh = ii
            break;
          }
        } else {
          miss = 0
        }
        pt.bitString[[sig.nodes[p]]][ii] = ind_t
      }
      pt.bitString[[sig.nodes[p]]] = pt.bitString[[sig.nodes[p]]][1:thresh]
      names(pt.bitString[[sig.nodes[p]]]) = perms[[sig.nodes[p]]][1:thresh]
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
    pt.byK[[k]] = pt.bitString[[which.min(bestInd)]]
  }
  
  return(pt.byK)
}


