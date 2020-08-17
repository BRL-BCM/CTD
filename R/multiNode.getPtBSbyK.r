#' Generate patient-specific bitstrings from adaptive network walk.
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss) associated with the adaptive network walk
#' made by the diffusion algorithm trying to find the variables in the encoded subset, given the background knowledge graph.
#' @param S - A character vector of node names describing the node subset to be encoded.
#' @param ranks - The list of node ranks calculated over all possible nodes in S, starting with each node in subset of interest.
#' @return pt.byK - a list of bitstrings, with the names of the list elements the node names of the encoded nodes
#' @export multiNode.getPtBSbyK
#' @examples
#' # Get patient bitstrings for the first 4 patients in the Miller et al 2015 dataset.
#' data("Miller2015")
#' data_mx = Miller2015[-grep("x - ", rownames(data_mx)),grep("IEM", colnames(Miller2015))]
#' data_mx = data_mx[,c(1:4)]
#' # Look at the top 15 metabolites for each patient. 
#' kmx=15
#' topMets_allpts = c()
#' for (pt in 1:ncol(data_mx)) { topMets_allpts = c(topMets_allpts, rownames(data_mx)[order(abs(data_mx[,pt]), decreasing=TRUE)[1:kmx]]) }
#' topMets_allpts = unique(topMets_allpts)
#' # Pre-compute node ranks for all metabolites in topMets_allpts for faster calculations.
#' ranks = list()
#' for (n in 1:length(topMets_allpts)) { ranks[n] = multiNode.getNodeRanksN(n, G, p1=0.9, thresholdDiff=0.01, adj_mat) }
#' names(ranks) = topMets_allpts
#' # Store each patient's bitstrings for each size k={1...kmx}.
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'   S = rownames(data_mx)[order(abs(data_mx[,pt]), decreasing=TRUE)[1:kmx]]
#'   ptBSbyK[[ptID]] = multiNode.getPtBSByK(S, ranks)
#' }
multiNode.getPtBSbyK = function(S, ranks) {
  ranks2 = ranks[which(names(ranks) %in% S)]
  pt.bitString = list()
  for (p in 1:length(ranks2)) {
    pt.bitString[[S[p]]] = as.numeric(ranks2[[S[p]]] %in% S)
    names(pt.bitString[[S[p]]]) = ranks2[[S[p]]]
    ind = which(pt.bitString[[S[p]]] == 1)
    pt.bitString[[S[p]]] = pt.bitString[[S[p]]][1:ind[length(ind)]]
  }
  # For each k, find the ranks that found at least k in S. Which found the first k soonest?
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
