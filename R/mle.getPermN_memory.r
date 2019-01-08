#' Generate the "adaptive walk" node permutations, starting from a given perturbed variable
#'
#' This function calculates the node permutation starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param n - The index (out of a vector of metabolite names) of the permutation you want to calculate.
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
mle.getPermN_memory = function(n, G) {
  all_nodes = names(G)
  print(sprintf("Calculating permutation %d of %d.", n, length(all_nodes)))
  current_node_set = NULL
  stopIterating=FALSE
  startNode = all_nodes[n]
  hits = startNode
  currentGraph = G
  while (stopIterating==FALSE) {
    current_node_set = c(current_node_set, startNode)
    sumHits = as.vector(matrix(0, nrow=length(G), ncol=1))
    names(sumHits) = names(G)
    for (hit in 1:length(hits)) {
      #For unseen nodes, clear probabilities and add probability (p0/#unseen nodes)
      for (t in 1:length(currentGraph)) {
        currentGraph[[t]] = 0; #set probabilities of all nodes to 0
      }
      #determine base p0 probability
      baseP = p0/(length(currentGraph)-length(current_node_set));
      for (t in 1:length(currentGraph)) {
        if (!(names(currentGraph[t]) %in% current_node_set)) {
          currentGraph[[t]] = baseP;  #set probabilities of unseen nodes to diffused p0 value, baseP
        } else {
          currentGraph[[t]] = 0;
        }
      }
      p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
      currentGraph = graph.diffuseP1(p1, hits[hit], currentGraph, currentGraph[current_node_set], 1, verbose=FALSE)
      p1_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
      if (abs(p1_event-1)>thresholdDiff) {
        extra.prob.to.diffuse = 1-p1_event
        currentGraph[names(current_node_set)] = 0
        currentGraph[!(names(currentGraph) %in% names(current_node_set))] = unlist(currentGraph[!(names(currentGraph) %in% names(current_node_set))]) + extra.prob.to.diffuse/sum(!(names(currentGraph) %in% names(current_node_set)))
      }
      sumHits = sumHits + unlist(currentGraph)
    }
    sumHits = sumHits/length(hits)
    #Set startNode to a node that is the max probability in the new currentGraph
    maxProb = names(which.max(sumHits))
    # Break ties: When there are ties, choose the first of the winners.
    startNode = names(currentGraph[maxProb[1]])
    if (length(c(startNode,current_node_set))>=(length(G))) {
      current_node_set = c(current_node_set, startNode)
      stopIterating = TRUE
    }
  }
  return(current_node_set)
}


