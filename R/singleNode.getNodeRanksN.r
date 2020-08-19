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
#' @param output_dir - If specified, a image sequence will generate in the output directory specified.
#' @param useLabels - If TRUE, node names will display next to their respective nodes in the network. If false
#'                    node names will not display. Only relevant if output_dir is specified. 
#' @param coords - The x and y coordinates for each node in the network, to remain static between images.
#' @return curr_nodeset - A character vector of node names in the order they were drawn by the probability diffusion algorithm.
#' @keywords probability diffusion
#' @keywords network walker
#' @export singleNode.getNodeRanksN
#' @examples
#' data("Miller2015")
#' data_mx = Miller2015[-c(1, grep("x - ", rownames(Miller2015))), grep("IEM", colnames(Miller2015))]
#' data_mx = apply(data_mx, c(1,2), as.numeric)
#' # Build an adjacency matrix for network G
#' adj_mat = matrix(0, nrow=nrow(data_mx), ncol=nrow(data_mx))
#' rows = sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' cols = sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' for (i in rows) {for (j in cols) { adj_mat[i, j] = rnorm(1, mean=0, sd=1)} }
#' colnames(adj_mat) = rownames(data_mx)
#' rownames(adj_mat) = rownames(data_mx)
#' G = vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' # Get node rankings for the first metabolite in network G. 
#' ranks = singleNode.getNodeRanksN(1, G, p1=0.9, thresholdDiff=0.01, adj_mat)
#' # Make a movie of the network walker
#' S = names(G)[sample(seq_len(length(G)), 3, replace=FALSE)]
#' ranks = singleNode.getNodeRanksN(which(names(G)==S[1]), G, p1=0.9, thresholdDiff=0.01, adj_mat, S, log2(length(G)), FALSE, getwd())
singleNode.getNodeRanksN = function(n, G, p1, thresholdDiff, adj_mat, S=NULL, num.misses=NULL, verbose=FALSE, output_dir="", useLabels=FALSE, coords=NULL) {
  p0 = 1 - p1
  if (is.null(S) && !is.null(num.misses)) {
    print("You must supply a subset of nodes as parameter S if you supply num.misses.")
    return(0)
  }
  if (is.null(S) && output_dir!="") {
    print("You must supple a subset of nodes as parameter S if you supply output_dir.")
    return(0)
  }
  if (verbose) { print(sprintf("Calculating node rankings %d of %d.", n, length(G))) }
  
  curr_nodeset = NULL
  stopIterating=FALSE
  startNode = names(G)[n]
  currGph = G
  numMisses = 0
  curr_nodeset = c(curr_nodeset, startNode)
  if (output_dir!="") { graph.takeNetWalkSnapShot(adj_mat, G, output_dir, p1, curr_nodeset, S, 
                                                  coords, imgNum=length(curr_nodeset), useLabels) }
  while (stopIterating==FALSE) {
    # Clear probabilities
    currGph[seq_len(length(currGph))] = 0 #set probabilities of all nodes to 0
    #determine base p0 probability
    baseP = p0/(length(currGph)-length(curr_nodeset))
    #set probabilities of unseen nodes to baseP
    currGph[!(names(currGph) %in% curr_nodeset)] = baseP
    # Sanity check. p0_event should add up to exactly p0 (global variable)
    p0_event = sum(unlist(currGph[!(names(currGph) %in% curr_nodeset)]))
    currGph = graph.diffuseP1(p1, startNode, currGph, curr_nodeset, thresholdDiff, adj_mat, verbose=FALSE)
    # Sanity check. p1_event should add up to exactly p1 (global variable)
    p1_event = sum(unlist(currGph[!(names(currGph) %in% curr_nodeset)]))
    if (abs(p1_event-1)>thresholdDiff) {
      extra.prob.to.diffuse = 1-p1_event
      currGph[names(curr_nodeset)] = 0
      currGph[!(names(currGph) %in% names(curr_nodeset))] = unlist(currGph[!(names(currGph) %in% names(curr_nodeset))]) + 
        extra.prob.to.diffuse/sum(!(names(currGph) %in% names(curr_nodeset)))
    }
    #Set startNode to a node that is the max probability in the new currGph
    maxProb = names(which.max(currGph))
    if (output_dir!="") { graph.takeNetWalkSnapShot(adj_mat, G, output_dir, p1, curr_nodeset, S, 
                                                    coords, imgNum=length(curr_nodeset), useLabels) }
    
    # Break ties: When there are ties, choose the first of the winners.
    startNode = names(currGph[maxProb[1]])
    if (!is.null(S)) {
      if (startNode %in% S) {
        numMisses = 0
      } else {
        numMisses = numMisses + 1
      }
      curr_nodeset = c(curr_nodeset, startNode)
      if (numMisses>num.misses || all(S %in% curr_nodeset)) {
        stopIterating = TRUE
      }
    } else {
      # Is S isn't specified, keep drawing until you've drawn all nodes.
      curr_nodeset = c(curr_nodeset, startNode)
      if (length(curr_nodeset)>=(length(G))) {
        stopIterating = TRUE
      }
    }
    if (output_dir!="") { graph.takeNetWalkSnapShot(adj_mat, G, output_dir, p1, curr_nodeset, S, 
                                                    coords, imgNum=length(curr_nodeset), useLabels) }
  }
  return(curr_nodeset)
}
