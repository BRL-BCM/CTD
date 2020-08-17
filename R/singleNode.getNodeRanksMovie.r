#' Capture the movement of the fixed, single-node walk of the diffusion probability method.
#'
#' Make a movie of the fixed, single-node walk the diffusion probability method makes in search of a given patient's perturbed variables.
#' @param S - The subset of variables, S, in a background graph, G.
#' @param ig - The igraph object associated with the network G.
#' @param output_filepath - The local directory at which you want still PNG images to be saved.
#' @param p1 - The probability that is preferentially distributed between network nodes by the 
#'             probability diffusion algorithm based solely on network connectivity. The remaining probability
#'             (i.e., "p0") is uniformally distributed between network nodes, regardless of connectivity.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges this amount (thresholdDiff)
#'                        or less between nodes, the algorithm returns up the call stack.
#' @param num.misses - The number of missteps allowed before a network walker truncates its walk.
#' @param zoomIn - Boolean. Delete nodes outside of node subset's order 1 neighborhood?. Default is FALSE.
#' @return ranksByStartNode - a list object of node rankings Each element is based on a different startNode.
#'         Images are also generated in the output_directory specified.
#' @importFrom igraph V E degree delete.vertices neighborhood get.adjacency graph.adjacency layout.fruchterman.reingold
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend title
#' @export singleNode.getNodeRanksMovie
#' @examples
#' # Read in any network via its adjacency matrix
#' adj_mat = matrix(1, nrow=100, ncol=100)
#' for (i in 1:100) {for (j in 1:100) { adj_mat[i, j] = rnorm(1, mean=0, sd=1)} }
#' colnames(adj_mat) = sprintf("Compound%d", 1:100)
#' ig = graph.adjacency(adj_mat, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = V(ig)$name
#' S = names(G)[sample(1:length(G), 3)]
#' singleNode.getNodeRanksMovie(S, ig, output_filepath = getwd(), p1=0.9, thresholdDiff=0.01)
singleNode.getNodeRanksMovie = function(S, ig, output_filepath, p1, thresholdDiff, num.misses=NULL, zoomIn=FALSE) {
  p0 = 1 - p1
  if (zoomIn) {
    # Experimental feature to "zoom in" on a node neighborhood in very large networks, to better visualize 
    # the node rankings instead of a hairball
    ig = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% names(unlist(neighborhood(ig, nodes=S, order=1))))])
    V(ig)$label.cex = 2
    V(ig)$label = rep("", length(V(ig)$name))
    tmp = get.adjacency(ig, attr="weight")
    tmp = abs(tmp)
    igraphTestG = graph.adjacency(tmp, mode="undirected", weighted=TRUE)
  }
  coords = layout.fruchterman.reingold(ig)
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name
  G = lapply(G, function(i) i[[1]]=0)
  degs = list(degree(ig))
  adj_mat = as.matrix(get.adjacency(ig, attr="weight"))

  # Do node rankings ahead of time for each possible startNode in subGraphS.
  V(ig)$color = rep("white", length(G))
  V(ig)$color[which(V(ig)$name %in% S)] = "green"
  V(ig)$label = sprintf("%s:%.2f", V(ig)$name, G)
  current_node_set = NULL
  png(sprintf("%s/diffusionP1Movie%d.png", output_filepath, length(current_node_set)), 500, 500)
  plot.igraph(ig, layout=coords, vertex.color=V(ig)$color, edge.arrow.size=0.01,
              vertex.label=V(ig)$label, edge.width=5*abs(E(ig)$weight),  mark.col="dark red",
              mark.groups = current_node_set)
  title(sprintf("{%s}", paste(as.numeric(current_node_set %in% S), collapse="")), cex.main=2)
  legend("topright", legend=c("Jackpot Nodes", "Drawn Nodes"), fill=c("green", "dark red"))
  dev.off()
  
  # Get all node names, so we know what all possible startNodes are.
  ranksByStartNode = list()
  bitStrings.pt = list()
  for (n in 1:length(S)) {
    print(sprintf("Calculating node rankings %d of %d...", n, length(S)))
    # Draw all nodes in graph
    current_node_set = NULL
    stopIterating=FALSE;
    startNode = S[n]
    hits = startNode
    numMisses = 0
    currentGraph = G
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

      V(ig)$label = sprintf("%s:%.2f", V(ig)$name, G)
      png(sprintf("%s/diffusionP1Movie%d_%d.png", output_filepath, n, length(current_node_set)), 500, 500)
      plot.igraph(ig, layout=coords, vertex.color=V(ig)$color,
                  vertex.label=V(ig)$label, edge.width=5*abs(E(ig)$weight), edge.arrow.size=0.01,
                  vertex.size=5+round(50*unlist(currentGraph), 0), mark.col="dark red",
                  mark.groups = current_node_set)
      title(sprintf("{%s}", paste(as.numeric(current_node_set %in% S), collapse="")), cex.main=2)
      legend("topright", legend=c("Jackpot Nodes", "Drawn Nodes"), fill=c("green", "dark red"))
      dev.off()
    }

    ranksByStartNode[[n]] = current_node_set
    bitStrings.pt[[n]] = as.numeric(current_node_set %in% S)
  }
  names(ranksByStartNode) = S

  return(ranksByStartNode)
}

