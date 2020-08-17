#' Generate the fixed , single-node diffusion node rankings, starting from a given perturbed variable.
#'
#' This function calculates the node rankings starting from a given perturbed variable in a subset of variables in the background knowledge graph.
#' @param n - The index (out of a vector of node names) of the node ranking you want to calculate.
#' @param G - A list of probabilities with list names being the node names of the background graph.
#' @param S - A character vector of node names in the subset you want the network walker to find.
#' @param num.misses - The number of "misses" the network walker will tolerate before switching to fixed length codes for remaining nodes to be found.
#' @param p1 - The probability that is preferentially distributed between network nodes by the 
#'             probability diffusion algorithm based solely on network connectivity. The remaining probability
#'             (i.e., "p0") is uniformally distributed between network nodes, regardless of connectivity.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges this amount (thresholdDiff)
#'                        or less between nodes, the algorithm returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for the network, G. 
#' @param verbose - If TRUE, print statements will execute as progress is made. Default is FALSE.
#' @return current_node_set - A character vector of node names in the order they were drawn by the probability diffusion algorithm.
#' @keywords probability diffusion
#' @keywords network walker
#' @export singleNode.getNodeRanksN
#' @examples
#' data("Miller2015")
#' data_mx = Miller2015[-c(1, grep("x - ", rownames(Miller2015))), grep("IEM", colnames(Miller2015))]
#' data_mx = apply(data_mx, c(1,2), as.numeric)
#' # Build an adjacency matrix for network G
#' adj_mat = matrix(0, nrow=nrow(data_mx), ncol=nrow(data_mx))
#' rows = sample(1:ncol(adj_mat), 0.1*ncol(adj_mat))
#' cols = sample(1:ncol(adj_mat), 0.1*ncol(adj_mat))
#' for (i in rows) {for (j in cols) { adj_mat[i, j] = rnorm(1, mean=0, sd=1)} }
#' colnames(adj_mat) = rownames(data_mx)
#' rownames(adj_mat) = rownames(data_mx)
#' G = vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' # Get node rankings for the first metabolite in network G. 
#' ranks = singleNode.getNodeRanksN(1, G, p1=0.9, thresholdDiff=0.01, adj_mat)
singleNode.getNodeRanksN = function(n, G, p1, thresholdDiff, adj_mat, S=NULL, num.misses=NULL, verbose=FALSE) {
  p0 = 1 - p1
  if (!is.null(num.misses)) {
    if (is.null(S)) {
      print("You must supply a subset of nodes as parameter S if you supply num.misses.")
      return(0)
    }
  }
  all_nodes = names(G)
  if (verbose) {
    print(sprintf("Calculating node rankings %d of %d.", n, length(all_nodes)))
  }
  current_node_set = NULL
  stopIterating=FALSE
  startNode = all_nodes[n]
  currentGraph = G
  numMisses = 0
  current_node_set = c(current_node_set, startNode)
  while (stopIterating==FALSE) {
    # Clear probabilities
    currentGraph[1:length(currentGraph)] = 0 #set probabilities of all nodes to 0
    #determine base p0 probability
    baseP = p0/(length(currentGraph)-length(current_node_set))
    #set probabilities of unseen nodes to baseP
    currentGraph[!(names(currentGraph) %in% current_node_set)] = baseP
    # Sanity check. p0_event should add up to exactly p0 (global variable)
    p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
    currentGraph = graph.diffuseP1(p1, startNode, currentGraph, current_node_set, thresholdDiff, adj_mat, verbose=FALSE)
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
    if (!is.null(num.misses)) {
      if (startNode %in% S) {
        numMisses = 0
      } else {
        numMisses = numMisses + 1
      }
      current_node_set = c(current_node_set, startNode)
      if (numMisses>num.misses || length(c(startNode,current_node_set))>=(length(G))) {
        stopIterating = TRUE
      }
    } else {
      # Keep drawing until you've drawn all nodes.
      current_node_set = c(current_node_set, startNode)
      if (length(current_node_set)>=(length(G))) {
        stopIterating = TRUE
      }
    }

  }
  return(current_node_set)
}
