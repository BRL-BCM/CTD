#' Generate the "adaptive walk" node rankings, starting from a given perturbed variable
#'
#' This function calculates the node rankings starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param S - A character vector of the node names for the subset of nodes you want to encode.
#' @param G - A list of probabilities with list names being the node names of the background graph.
#' @param p1 - The probability that is preferentially distributed between network nodes by the 
#'             probability diffusion algorithm based solely on network connectivity. The remaining probability
#'             (i.e., "p0") is uniformally distributed between network nodes, regardless of connectivity.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges this amount (thresholdDiff)
#'                        or less between nodes, the algorithm returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for the network, G. 
#' @param num.misses - The number of "misses" the network walker will tolerate before switching to fixed length codes for remaining nodes to be found.
#' @param verbose - If TRUE, print statements will execute as progress is made. Default is FALSE.
#' @return ranks - A list of character vectors of node names in the order they were drawn by the
#'                 probability diffusion algorithm, from each starting node in S.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export multiNode.getNodeRanks
#' @examples
#' # Read in any network via its adjacency matrix
#' adj_mat = matrix(1, nrow=100, ncol=100)
#' for (i in seq_len(100)) {for (j in seq_len(100)) {adj_mat[i, j] = rnorm(1, mean=0, sd=1)}}
#' colnames(adj_mat) = sprintf("Metabolite%d", seq_len(100))
#' rownames(adj_mat) = colnames(adj_mat)
#' G = vector(mode="list", length=ncol(adj_mat))
#' names(G) = colnames(adj_mat)
#' S = names(G)[seq_len(3)]
#' ranks = multiNode.getNodeRanks(S, G, p1=0.9, thresholdDiff=0.01, adj_mat)
multiNode.getNodeRanks = function(S, G, p1, thresholdDiff, adj_mat, num.misses=NULL, verbose=FALSE) {
  p0 = 1-p1
  if (is.null(num.misses)) {
    num.misses = log2(length(G))
  }
  ranks = list()
  for (n in seq_len(length(S))) {
    if (verbose) {
      print(sprintf("Calculating node rankings %d of %d.", n, length(S)))
    }
    current_node_set = NULL
    stopIterating=FALSE
    startNode = S[n]
    hits = startNode
    numMisses = 0
    current_node_set = c(current_node_set, startNode)
    while (stopIterating==FALSE) {
      sumHits = as.vector(matrix(0, nrow=length(G), ncol=1))
      names(sumHits) = names(G)
      for (hit in seq_len(length(hits))) {
        #Clear probabilities.
        G[seq_len(length(G))] = 0
        #determine base p0 probability
        baseP = p0/(length(G)-length(current_node_set))
        #set probabilities of unseen nodes to diffused p0 value, baseP
        G[which(!(names(G) %in% current_node_set))] = baseP
        # Sanity check: p0_event should be equal to p0 global variable.
        p0_event = sum(unlist(G[!(names(G) %in% current_node_set)]))
        # Diffuse p1.
        G = graph.diffuseP1(p1, hits[hit], G, current_node_set, thresholdDiff, adj_mat, verbose=FALSE)
        # Sanity check: p1_event should be within 'thresholdDiff' (global variable) of 1.
        p1_event = sum(unlist(G[!(names(G) %in% current_node_set)]))
        if (abs(p1_event-1)>thresholdDiff) {
          extra.prob.to.diffuse = 1-p1_event
          G[names(current_node_set)] = 0
          unseen_nodes = which(!(names(G) %in% names(current_node_set)))
          G[unseen_nodes] = unlist(G[unseen_nodes]) + extra.prob.to.diffuse/length(unseen_nodes)
        }
        sumHits = sumHits + unlist(G)
      }
      sumHits = sumHits/length(hits)
      #Set startNode to a node that is the max probability in the new G
      maxProb = names(which.max(sumHits[-which(names(sumHits) %in% current_node_set)]))
      # Break ties: When there are ties, choose the first of the winners.
      startNode = names(G[maxProb[1]])
      if (startNode %in% S) {
        numMisses = 0
        #print(sprintf("Drew %s on draw %d, a hit!", startNode, length(current_node_set)+1))
        hits = c(hits, startNode)
      } else {
        numMisses = numMisses + 1
      }

      current_node_set = c(current_node_set, startNode)
      if (all(S %in% current_node_set) || numMisses>num.misses) {
        stopIterating = TRUE
      }

    }
    ranks[[n]] = current_node_set
  }
  names(ranks) = S
  return(ranks)
}


