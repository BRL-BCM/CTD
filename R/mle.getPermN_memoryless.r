#' Generate the "adaptive walk" node permutations, starting from a given perturbed variable
#'
#' This function calculates the node permutation starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param n - The index (out of a vector of metabolite names) of the permutation you want to calculate.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export mle.getPermN
#' @examples
#' # Read in any network via its adjacency matrix
#' tmp = matrix(1, nrow=100, ncol=100)
#' for (i in 1:100) {
#'   for (j in 1:100) {
#'     tmp[i, j] = rnorm(1, mean=0, sd=1)
#'   }
#' }
#' colnames(tmp) = sprintf("MolPheno%d", 1:100)
#' ig = graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))  # Must have this declared as a GLOBAL variable!!!!!
#' # Set other tuning parameters
#' p0=0.1  # 10% of probability distributed uniformly
#' p1=0.9  # 90% of probability diffused based on edge weights in networks
#' thresholdDiff=0.01
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = V(ig)$name
#' # Get node permutations for graph
#' perms = list()
#' for (n in 1:length(G)) {
#'   print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
#'   perms[[n]] = mle.getPermN(n, G)
#' }
#' names(perms) = names(G)
mle.getPermN = function(n, G) {
  all_nodes = names(G)
  print(sprintf("Calculating permutation %d of %d.", n, length(all_nodes)))
  current_node_set = NULL
  stopIterating=FALSE
  startNode = all_nodes[n]
  currentGraph = G
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
    if (length(c(startNode,current_node_set))>=(length(G))) {
      current_node_set = c(current_node_set, startNode)
      stopIterating = TRUE
    }
  }

  return(current_node_set)
}


