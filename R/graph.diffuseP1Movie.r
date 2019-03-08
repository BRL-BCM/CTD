#' Make a movie of the diffusion of probability, P1, from a starting node.
#'
#' Recursively diffuse probability from a starting node based on the connectivity of the background knowledge graph, representing the likelihood that a variable will be
#'         most influenced by a perturbation in the starting node.
#' @param p1 - The probability being dispersed from the starting node, startNode.
#' @param startNode - The first variable drawn in the adaptive permutation node sequence, from which p1 gets dispersed.
#' @param G - A list of probabilities, with names of the list being the node names in the background knowledge graph.
#' @param visitedNodes - A character vector of node names, storing the history of previous draws in the permutation sequence.
#' @param graphNumber - If testing against multiple background knowledge graphs, this is the index associated with the adjacency matrix that codes for G. Default value is 1.
#' @return G - A list of returned probabilities after the diffusion of probability has truncated, with names of the list being the node names in the background knowledge graph.
#' @export graph.diffuseP1Movie
#' @keywords generative methods
#' @keywords diffusion event
#' @keywords network walker
#' @examples
#' # 7 node example graph from Figure 4 (Thistlethwaite, Elsea, Milosavljevic, 2019)
#' adj_mat = rbind(c(0,2,1,0,0,0,0), # A
#'                 c(2,0,1,0,0,0,0), # B
#'                 c(1,0,0,1,0,0,0), # C
#'                 c(0,0,1,0,2,0,0), # D
#'                 c(0,0,0,2,0,2,1), # E
#'                 c(0,0,0,1,2,0,1), # F
#'                 c(0,0,0,0,1,1,0)  # G
#'                 )
#' rownames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G")
#' colnames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G")
#' ig = graph.adjacency(as.matrix(adj_mat), mode="undirected", weighted=TRUE)
#' G=vector(mode="list", length=7)
#' G[1:length(G)] = 0
#' names(G) = c("A", "B", "C", "D", "E", "F", "G")
#' startNode = "A"
#' visitedNodes = startNode
#' # Diffuse 100% of probability from startNode "A"
#' p1 = 1.0
#' # Probability diffusion truncates at
#' thresholdDiff=0.01
#' coords = layout.fruchterman.reingold(ig)
#' V(ig)$x = coords[,1]
#' V(ig)$y = coords[,2]
#' # Global variable imgNum
#' imgNum=1
#' G_new = graph.diffuseP1Movie(p1, startNode, G, visitedNodes, ig, 1, getwd())
graph.diffuseP1Movie = function(p1, startNode, G, visitedNodes, ig, recursion_level=1, output_dir=getwd()) {
  print(sprintf("%sprob. to diffuse:%f startNode: %s, visitedNodes: %s",
                paste(rep("   ", length(visitedNodes) - 1), collapse = ""),
                p1, startNode, toString(visitedNodes)))
  print(sprintf("ImgNum=%d", imgNum))

  V(ig)$color = rep("blue", length(G))
  V(ig)$color[which(V(ig)$name %in% visitedNodes)] = "red"
  V(ig)$label = sprintf("%s:%.2f", V(ig)$name, G)
  png(sprintf("%s/diffusionEventMovie%d.png", output_dir, .GlobalEnv$imgNum), 500, 500)
  plot.igraph(ig, layout=cbind(V(ig)$x, V(ig)$y), vertex.color=V(ig)$color,
              vertex.label=V(ig)$label, vertex.label.dist = 3, edge.width=5*abs(E(ig)$weight),
              mark.col="black", mark.border = "black", mark.groups = startNode)
  title(sprintf("Diffuse %.2f from %s at recursion level %d.", p1, startNode, recursion_level), cex.main=1)
  legend("bottomright", legend=c("Visited", "Unvisited"), fill=c("red", "blue"))
  dev.off()
  .GlobalEnv$imgNum=.GlobalEnv$imgNum+1

  adj_mat = as.matrix(get.adjacency(ig, attr = "weight"))
  startNodeNeighbors = names(which(abs(adj_mat[, startNode]) >  0))
  startNodeUnvisitedNeighbors = startNodeNeighbors[!(startNodeNeighbors %in% visitedNodes)]
  vN = visitedNodes[which(visitedNodes != startNode)]
  extendedConnections = NULL
  if (length(vN) > 0 && length(startNodeUnvisitedNeighbors) == 0) {
    adj_matAfter = adj_mat[-which(rownames(adj_mat) %in% vN), -which(colnames(adj_mat) %in% vN)]
    connections = adj_mat[startNode, ]
    connectionsYes = connections[which(abs(connections) > 0)]
    connectionsNo = connections[intersect(which(connections == 0), which(!(names(connections) %in% c(startNode, vN))))]
    if (length(connectionsNo) > 0) {
      for (n1 in 1:length(connectionsNo)) {
        if (length(connectionsYes) > 0) {
          for (n2 in 1:length(connectionsYes)) {
            if (abs(adj_mat[names(connectionsYes[n2]), names(connectionsNo[n1])]) > 0) {
              connectionsNo[n1] = adj_mat[names(connectionsYes[n2]), names(connectionsNo[n1])]
              extendedConnections = c(extendedConnections, connectionsNo[n1])
            }
          }
        }
      }
    }
    if (length(extendedConnections) > 0) {
      adj_matAfter[startNode, names(extendedConnections)] = extendedConnections
      adj_matAfter[names(extendedConnections), startNode] = extendedConnections
    }
    adj_mat = adj_matAfter
  }
  startNodeNeighbors = names(which(abs(adj_mat[, startNode]) > 0))
  startNodeUnvisitedNeighbors = startNodeNeighbors[!(startNodeNeighbors %in% visitedNodes)]
  if (length(startNodeUnvisitedNeighbors) > 0 || length(extendedConnections) > 0) {
    weighted_edges.sum = sum(abs(adj_mat[which(rownames(adj_mat) %in% startNodeUnvisitedNeighbors), startNode]))
    print(sprintf("%sWeighted_edges.sum=%f", paste(rep("   ", length(visitedNodes) - 1), collapse = ""), weighted_edges.sum))
    z = 1
    while (z <= length(startNodeUnvisitedNeighbors)) {
      inherited.probability = p1 * abs(adj_mat[startNodeUnvisitedNeighbors[z], startNode])/weighted_edges.sum
      G[[startNodeUnvisitedNeighbors[z]]] = G[[startNodeUnvisitedNeighbors[z]]] + inherited.probability

      print(sprintf("%sadded diffused  probability %f to child #%d: %s", paste(rep("   ", length(visitedNodes) - 1), collapse = ""), inherited.probability, z, startNodeUnvisitedNeighbors[z]))

      V(ig)$color = rep("blue", length(G))
      V(ig)$color[which(V(ig)$name %in% visitedNodes)] = "red"
      V(ig)$label = sprintf("%s:%.2f", V(ig)$name, G)
      #tmp = rep(0.01, length(E(ig)$weight))
      #tmp[get.edge.ids(ig, vp=c(startNode, startNodeUnvisitedNeighbors[z]))] = 1
      lbls = rep("", length(E(ig)$weight))
      lbls[get.edge.ids(ig, vp=c(startNode, startNodeUnvisitedNeighbors[z]))] = round(inherited.probability, 3)
      png(sprintf("%s/diffusionEventMovie%d.png", output_dir, .GlobalEnv$imgNum), 500, 500)
      plot.igraph(ig, layout=cbind(V(ig)$x, V(ig)$y), vertex.color=V(ig)$color,
                  vertex.label=V(ig)$label, vertex.label.dist = 3, edge.label = lbls, edge.width=5*abs(E(ig)$weight),
                  mark.col="black", mark.border = "black", mark.groups = startNode)
      title(sprintf("Diffuse %.2f from %s at recursion level %d.", p1, startNode, recursion_level), cex.main=1)
      legend("bottomright", legend=c("Visited", "Unvisited"), fill=c("red", "blue"))
      dev.off()
      .GlobalEnv$imgNum = .GlobalEnv$imgNum + 1

      nNeighbors = names(which(abs(adj_mat[, startNodeUnvisitedNeighbors[z]]) >  0))
      if (length(nNeighbors) > 0 && inherited.probability/2 > thresholdDiff && ((length(visitedNodes) + 1) < length(G))) {
        G[[startNodeUnvisitedNeighbors[z]]] = G[[startNodeUnvisitedNeighbors[z]]] - inherited.probability/2

        print(sprintf("%ssubtracted %f from startNode's neighbor #%d: %s and sent to its own neighbors.",
                        paste(rep("   ", length(visitedNodes) - 1), collapse = ""), inherited.probability/2,
                        z, startNodeUnvisitedNeighbors[z]))

        V(ig)$color = rep("blue", length(G))
        V(ig)$color[which(V(ig)$name %in% visitedNodes)] = "red"
        V(ig)$label = sprintf("%s:%.2f", V(ig)$name, G)
        #tmp = rep(0.01, length(E(ig)$weight))
        #tmp[get.edge.ids(ig, vp=c(startNode, startNodeUnvisitedNeighbors[z]))] = 1
        lbls = rep("", length(E(ig)$weight))
        lbls[get.edge.ids(ig, vp=c(startNode, startNodeUnvisitedNeighbors[z]))] = round(inherited.probability/2, 3)
        png(sprintf("%s/diffusionEventMovie%d.png", output_dir, .GlobalEnv$imgNum), 500, 500)
        plot.igraph(ig, layout=cbind(V(ig)$x, V(ig)$y), vertex.color=V(ig)$color,
                    vertex.label=V(ig)$label, vertex.label.dist = 3, edge.label = lbls, edge.width=5*abs(E(ig)$weight),
                    mark.col="black", mark.border = "black", mark.groups = startNode)
        title(sprintf("Diffuse %.2f from %s at recursion level %d.", p1, startNode, recursion_level), cex.main=1)
        legend("bottomright", legend=c("Visited", "Unvisited"), fill=c("red", "blue"))
        dev.off()
        .GlobalEnv$imgNum = .GlobalEnv$imgNum + 1

        G = graph.diffuseP1Movie(inherited.probability/2, startNodeUnvisitedNeighbors[z], G,
                                 c(visitedNodes, startNodeUnvisitedNeighbors[z]), ig, recursion_level+1, output_dir)
      }
      z = z + 1
    }
    V(ig)$color = rep("blue", length(G))
    V(ig)$color[which(V(ig)$name %in% visitedNodes)] = "red"
    V(ig)$label = sprintf("%s:%.2f", V(ig)$name, G)
    png(sprintf("%s/diffusionEventMovie%d.png", output_dir, .GlobalEnv$imgNum), 500, 500)
    plot.igraph(ig, layout=cbind(V(ig)$x, V(ig)$y), vertex.color=V(ig)$color,
                vertex.label=V(ig)$label, vertex.label.dist = 3, edge.width=5*abs(E(ig)$weight),
                mark.col="black", mark.border = "black", mark.groups = startNode)
    title(sprintf("Diffuse %.2f from %s at recursion level %d.", p1, startNode, recursion_level), cex.main=1)
    legend("bottomright", legend=c("Visited", "Unvisited"), fill=c("red", "blue"))
    dev.off()
    .GlobalEnv$imgNum=.GlobalEnv$imgNum+1
  } else {
    print("startNode is stranded with its visited neighbors, or is a singleton. Diffuse p1 uniformly amongst all unvisited nodes")
    t = 1
    while (t <= length(G)) {
      if (!(names(G[t]) %in% visitedNodes)) {
        G[[t]] = G[[t]] + p1/(length(G) - length(visitedNodes))
      }
      t = t + 1
    }
  }
  return(G)
}
