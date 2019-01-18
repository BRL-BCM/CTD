#' Generate the "adaptive walk" node permutations, starting from a given perturbed variable
#'
#' This function calculates the node permutation starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param S - A character vector of the node names for the subset of nodes you want to encode.
#' @return current_node_set - A character vector of node names in the order they were drawn by the probability diffusion algorithm.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export mle.getPermN_memory
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Get node permutations for graph
#' perms = list()
#' for (n in 1:length(G)) {
#'   print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
#'   perms[[n]] = mle.getPermN_memory(n, G)
#' }
#' names(perms) = names(G)
mle.getPermN_memory = function(S, G) {
  thresholdDrawT = log2(length(G))
  perms = list()
  for (n in 1:length(S)) {
    print(sprintf("Calculating permutation %d of %d.", n, length(S)))
    current_node_set = NULL
    stopIterating=FALSE
    startNode = S[n]
    hits = startNode
    numMisses = 0
    while (stopIterating==FALSE) {
      current_node_set = c(current_node_set, startNode)
      sumHits = as.vector(matrix(0, nrow=length(G), ncol=1))
      names(sumHits) = names(G)
      for (hit in 1:length(hits)) {
        #Clear probabilities.
        G[1:length(G)] = 0
        #determine base p0 probability
        baseP = p0/(length(G)-length(current_node_set))
        #set probabilities of unseen nodes to diffused p0 value, baseP
        G[which(!(names(G) %in% current_node_set))] = baseP
        # Sanity check: p0_event should be equal to p0 global variable.
        p0_event = sum(unlist(G[!(names(G) %in% current_node_set)]))
        # Diffuse p1.
        G = graph.diffuseP1(p1, hits[hit], G, current_node_set, 1, verbose=FALSE)
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
        print(sprintf("Drew %s on draw %d, a hit!", startNode, length(current_node_set)+1))
        hits = c(hits, startNode)
      } else {
        numMisses = numMisses + 1
      }

      if (all(S %in% current_node_set) || numMisses>thresholdDrawT) {
        current_node_set = c(current_node_set, startNode)
        stopIterating = TRUE
      }
    }
    perms[[n]] = current_node_set
  }
  names(perms) = S
  return(perms)
}


