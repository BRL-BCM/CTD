#' Apply the Kraft-McMillian Inequality using a specific encoding algorithm.
#'
#' A power analysis of the encoding algorithm using to encode subsets of S in G.
#' @param G - A character vector of all node names in the background knowledge graph.
#' @param k - The size of the node name subsets of G.
#' @return IA - a list of bitlengths associated with all outcomes in the N choose K outcome space, with the names of the list elements the node names of the encoded nodes
#' @export mle.kraftMcMillian
#' @example
#' G = list(A=0, B=0, C=0, D=0, E=0, F=0, G=0)
#' names(G) = tolower(names(G))
#' adj_mat = rbind(c(0,2,1,0,0,0,0), #A's neighbors
#'                 c(2,0,1,0,0,0,0), #B's neighbors
#'                 c(1,1,0,1,0,0,0), #C's neighbors
#'                 c(0,0,1,0,2,1,0), #D's neighbors
#'                 c(0,0,0,2,0,2,1), #E's neighbors
#'                 c(0,0,0,1,2,0,1), #F's neighbors
#'                 c(0,0,0,0,1,1,0)  #G's neighbors
#'                 )
#' rownames(adj_mat) = names(G)
#' colnames(adj_mat) = names(G)
#' adjacency_matrix = list(adj_mat)
#' IA = mle.kraftMcMillian(G, 2)
#' # Power to find effects is
#' sum(2^-unlist(IA))
mle.kraftMcMillian = function(G, k, memory=FALSE) {
  IA = list()
  # Get all subsets of size k in graph G
  subsets.k = combn(tolower(names(G)),k)
  for (ss in 1:ncol(subsets.k)) {
    S = subsets.k[,ss]
    perms = list()
    for (i in 1:length(S)) {
      ind = which(names(G)==S[i])
      if (memory) {
        perms[[S[i]]] = mle.getPermN_memory(ind, G)
      } else {
        perms[[S[i]]] = mle.getPermN_memoryless(ind, G)
      }
    }
    pt.bitString = list()
    for (p in 1:length(S)) {
      pt.bitString[[S[p]]] = as.numeric(perms[[S[p]]] %in% S)
      names(pt.bitString[[S[p]]]) = perms[[S[p]]]
      ind = which(pt.bitString[[S[p]]] == 1)
      pt.bitString[[S[p]]] = pt.bitString[[S[p]]][1:ind[length(ind)]]
    }
    # Which found the nodes in S soonest
    bestInd = vector("numeric", length(pt.bitString))
    for (p in 1:length(pt.bitString)) {
      bestInd[p] = sum(which(pt.bitString[[p]] == 1))
    }
    optBS = pt.bitString[[which.min(bestInd)]]
    mets.k = names(optBS)[which(optBS==1)]
    # Encoding length of subset S in BITS
    IA[[ss]] = log2(length(G)) + length(optBS)-1
    names(IA[[ss]]) = paste(S, collapse="")
  }

  return(IA)
}


