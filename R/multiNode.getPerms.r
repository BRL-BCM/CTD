#' Generate the "adaptive walk" node permutations, starting from a given perturbed variable
#'
#' This function calculates the node permutation starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param S - A character vector of the node names for the subset of nodes you want to encode.
#' @return current_node_set - A character vector of node names in the order they were drawn by the probability diffusion algorithm.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export multiNode.getPerms
#' @examples
#' # Read in any network via its adjacency matrix
#' tmp = matrix(1, nrow=100, ncol=100)
#' for (i in 1:100) {
#'   for (j in 1:100) {
#'     tmp[i, j] = rnorm(1, mean=0, sd=1)
#'   }
#' }
#' colnames(tmp) = sprintf("Compound%d", 1:100)
#' ig = graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' # Get node permutations for graph
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = V(ig)$name
#' S = names(G)[1:3]
#' perms = multiNode.getPerms(S, G)
multiNode.getPerms = function(S, G, num.misses=NULL) {
  if (is.null(num.misses)) {
    thresholdDrawT = log2(length(G))
  } else {
    thresholdDrawT = num.misses
  }
  perms = list()
  for (n in 1:length(S)) {
    print(sprintf("Calculating permutation %d of %d.", n, length(S)))
    current_node_set = NULL
    stopIterating=FALSE
    startNode = S[n]
    hits = startNode
    numMisses = 0
    current_node_set = c(current_node_set, startNode)
    while (stopIterating==FALSE) {
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
        #print(sprintf("Drew %s on draw %d, a hit!", startNode, length(current_node_set)+1))
        hits = c(hits, startNode)
      } else {
        numMisses = numMisses + 1
      }

      current_node_set = c(current_node_set, startNode)
      if (all(S %in% current_node_set) || numMisses>thresholdDrawT) {
        stopIterating = TRUE
      }

    }
    perms[[n]] = current_node_set
  }
  names(perms) = S
  return(perms)
}


