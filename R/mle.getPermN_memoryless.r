#' Generate the "adaptive walk" node permutations, starting from a given perturbed variable
#'
#' This function calculates the node permutation starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param n - The index (out of a vector of node names) of the permutation you want to calculate.
#' @param G - A list of probabilities with list names being the node names of the background graph.
#' @return current_node_set - A character vector of node names in the order they were drawn by the probability diffusion algorithm.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export mle.getPermN_memoryless
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Get node permutations for graph
#' perms = list()
#' for (n in 1:length(G)) {
#'   print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
#'   perms[[n]] = mle.getPermN(n, G)
#' }
#' names(perms) = names(G)
mle.getPermN_memoryless = function(n, G, S=NULL, misses.thresh=NULL) {
  if (!is.null(misses.thresh)) {
    if (is.null(S)) {
      print("You must supply a subset of nodes as parameter S if you supply a misses.thresh.")
      return(0)
    }
  }
  all_nodes = names(G)
  print(sprintf("Calculating permutation %d of %d.", n, length(all_nodes)))
  current_node_set = NULL
  stopIterating=FALSE
  startNode = all_nodes[n]
  currentGraph = G
  numMisses = 0
  while (stopIterating==FALSE) {
    current_node_set = c(current_node_set, startNode)
    # Clear probabilities
    currentGraph[1:length(currentGraph)] = 0 #set probabilities of all nodes to 0
    #determine base p0 probability
    baseP = p0/(length(currentGraph)-length(current_node_set))
    #set probabilities of unseen nodes to baseP
    currentGraph[!(names(currentGraph) %in% current_node_set)] = baseP
    # Sanity check. p0_event should add up to exactly p0 (global variable)
    p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
    currentGraph = graph.diffuseP1(p1, startNode, currentGraph, current_node_set, 1, verbose=FALSE)
    # Sanity check. p1_event should add up to exactly p1 (global variable)
    p1_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
    if (abs(p1_event-1)>thresholdDiff) {
      extra.prob.to.diffuse = 1-p1_event
      currentGraph[names(current_node_set)] = 0
      currentGraph[!(names(currentGraph) %in% names(current_node_set))] = unlist(currentGraph[!(names(currentGraph) %in% names(current_node_set))]) + extra.prob.to.diffuse/sum(!(names(currentGraph) %in% names(current_node_set)))
    }
    #Set startNode to a node that is the max probability in the new currentGraph
    maxProb = names(which.max(currentGraph))
    # Break ties: When there are ties, choose the first of the winners.
    startNode = names(currentGraph[maxProb[1]])
    if (!is.null(misses.thresh)) {
      if (startNode %in% S) {
        numMisses = 0
      } else {
        numMisses = numMisses + 1
      }
      if (numMisses>misses.thresh || length(c(startNode,current_node_set))>=(length(G))) {
        current_node_set = c(current_node_set, startNode)
        stopIterating = TRUE
      }
    } else {
      if (length(c(startNode,current_node_set))>=(length(G))) {
        current_node_set = c(current_node_set, startNode)
        stopIterating = TRUE
      }
    }

  }
  return(current_node_set)
}


mle.getPermN_trunLess = function(S, G, miss.thresh=50) {
  perms = list()
  for (n in 1:length(S)) {
    print(sprintf("Calculating permutation %d of %d.", n, length(S)))
    current_node_set = NULL
    stopIterating=FALSE
    startNode = S[n]
    currentGraph = G
    numMisses = 0
    while (stopIterating==FALSE) {
      current_node_set = c(current_node_set, startNode)
      # Clear probabilities
      currentGraph[1:length(currentGraph)] = 0 #set probabilities of all nodes to 0
      #determine base p0 probability
      baseP = p0/(length(currentGraph)-length(current_node_set))
      #set probabilities of unseen nodes to baseP
      currentGraph[!(names(currentGraph[t]) %in% current_node_set)] = baseP
      # Sanity check. p0_event should add up to exactly p0 (global variable)
      p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
      currentGraph = graph.diffuseP1(p1, startNode, currentGraph, currentGraph[current_node_set], 1, verbose=FALSE)
      # Sanity check. p1_event should add up to exactly p1 (global variable)
      p1_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
      if (abs(p1_event-1)>thresholdDiff) {
        extra.prob.to.diffuse = 1-p1_event
        currentGraph[names(current_node_set)] = 0
        ind = which(!(names(currentGraph) %in% names(current_node_set)))
        currentGraph[ind] = unlist(currentGraph[ind]) + extra.prob.to.diffuse/length(ind)
      }
      #Set startNode to a node that is the max probability in the new currentGraph
      maxProb = names(which.max(currentGraph))
      # Break ties: When there are ties, choose the first of the winners.
      startNode = names(currentGraph[maxProb[1]])
      if (startNode %in% S) {
        numMisses = 0
        print(sprintf("Drew %s on draw %d, a hit!", startNode, length(current_node_set)+1))
      } else {
        numMisses = numMisses + 1
      }

      if (all(S %in% current_node_set) || numMisses>miss.thresh) {
        current_node_set = c(current_node_set, startNode)
        stopIterating = TRUE
      }
    }
    perms[[n]] = current_node_set
  }
  names(perms) = S
  return(perms)
}


