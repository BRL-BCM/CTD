#' Generate the "adaptive walk" node permutations, starting from a given perturbed variable
#'
#' This function calculates the node permutation starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param n - The index (out of a vector of metabolite names) of the permutation you want to calculate.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export mle.getPermN
#' @examples
#' mle.getPermN(n)
mle.getPermN = function(n) {
  all_nodes = names(igraphObjectG)
  print(sprintf("Calculating permutation %d of %d.", n, length(all_nodes)))
  current_node_set = NULL
  stopIterating=FALSE
  startNode = all_nodes[n]
  currentGraph = igraphObjectG
  while (stopIterating==FALSE) {
    current_node_set = c(current_node_set, startNode)
    #STEP 4: Diffuse p0 and p1 to connected nodes from current draw, startNode node
    #For unseen nodes, clear probabilities and add probability (p0/#unseen nodes)
    for (t in 1:length(currentGraph)) {
      currentGraph[[t]] = 0 #set probabilities of all nodes to 0
    }
    #determine base p0 probability
    baseP = p0/(length(currentGraph)-length(current_node_set))
    for (t in 1:length(currentGraph)) {
      if (!(names(currentGraph[t]) %in% current_node_set)) {
        currentGraph[[t]] = baseP  #set probabilities of unseen nodes to diffused p0 value, baseP
      } else {
        currentGraph[[t]] = 0
      }
    }
    p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
    currentGraph = graph.diffuseP1(p1, startNode, currentGraph, currentGraph[current_node_set], 1, verbose=FALSE)
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
    if (length(c(startNode,current_node_set))>=(length(V(igraphTestG)$name))) {
      current_node_set = c(current_node_set, startNode)
      stopIterating = TRUE
    }
  }

  return(current_node_set)
}


