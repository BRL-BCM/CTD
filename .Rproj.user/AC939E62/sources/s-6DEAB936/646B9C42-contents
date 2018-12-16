#' Diffuse Probability P1 from a starting node.
#'
#' Recursively diffuse probability from a starting node based on the connectivity of the background knowledge graph, representing the likelihood that a variable will be
#'         most influenced by a perturbation in the starting node.
#' @param p1 - The probability being dispersed from the starting node, startNode.
#' @param startNode - The first variable drawn in the adaptive permutation node sequence, from which p1 gets dispersed.
#' @param G - The igraph object associated with the background knowledge graph.
#' @param visitedNodes - The history of previous draws in the permutation sequence.
#' @param graphNumber - If testing against multiple background knowledge graphs, this is the index associated with the adjacency matrix that codes for G. Default value is 1.
#' @param verbose - If debugging or tracking a diffusion event, verbose=TRUE will activate print statements. Default is FALSE.
#' @export graph.diffuseP1
#' @keywords generative methods
#' @keywords diffusion event
#' @keywords network walker
#' @examples
#' # Read in any network via its adjacency matrix
#' tmp = as.matrix(read.table("adjacency_matrix.txt", sep="\t", header=TRUE))
#' colnames(tmp) = rownames(tmp)
#' ig = graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))  # Must have this declared as a GLOBAL variable!!!!!
#'
#' p0=0.1  # 10% of probability distributed uniformly
#' p1=0.9  # 90% of probability diffused based on edge weights in networks
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = names(V(ig)$name)
#' startNode = names(G)[1]
#' visitedNodes = NULL
#' probs_afterCurrDraw = graph.diffuseP1(p1, startNode, G, visitedNodes, 1)
graph.diffuseP1 = function(p1, startNode, G, visitedNodes, graphNumber=1, verbose=FALSE) {
  if (verbose==TRUE) {
    print(sprintf("%sprob. to diffuse:%f startNode: %s, visitedNodes: %s", paste(rep("   ", length(visitedNodes)-1), collapse=""),
                  p1, startNode, toString(names(visitedNodes))));
  }
  adj_mat = adjacency_matrix[[graphNumber]];
  startNodeNeighbors = G[names(which(abs(adj_mat[,startNode]) > 0))]
  startNodeUnvisitedNeighbors = startNodeNeighbors[!(names(startNodeNeighbors) %in% names(visitedNodes))]

  # Diffuse only if startNode has unvisited neighbors.
  # Alternatively, startNode could have no unvisited neighbors, but before we decide to disperse p1 uniformly,
  # we need to ask if the visited neighbors of startNode have unvisited neighbors ("extended connections"),
  # or is startNode stranded alone with its visited neighbors?
  # Remove visited neighbors of startNode from adj_mat
  vN = visitedNodes[which(names(visitedNodes)!=startNode)]
  extendedConnections = NULL
  if (length(vN)>0 && length(startNodeUnvisitedNeighbors)==0) {
    adj_matAfter = adj_mat[-which(rownames(adj_mat) %in% names(vN)),
                            -which(colnames(adj_mat) %in% names(vN))]

    # Get extended unvisited neighbors of startNode from adj_mat
    # Get row in adj_mat for startNode
    connections = adj_mat[startNode,]
    connectionsYes = connections[which(abs(connections)>0)] #This should be 0 always??
    connectionsNo = connections[intersect(which(connections==0),
                                           which(!(names(connections) %in% c(startNode, names(vN)))))]
    if (length(connectionsNo)>0) {
      for (n1 in 1:length(connectionsNo)) {
        # For any 0 in connections, is it an unvisited node that is 1 in any of connection's 1's?
        if (length(connectionsYes)>0) {
          for (n2 in 1:length(connectionsYes)) {
            if (abs(adj_mat[names(connectionsYes[n2]),names(connectionsNo[n1])])>0) {
              connectionsNo[n1] = adj_mat[names(connectionsYes[n2]),names(connectionsNo[n1])]
              extendedConnections = c(extendedConnections, connectionsNo[n1])
              #print(sprintf("%s extended connection of %s via %s", names(connectionsNo[n1]), startNode, names(connectionsYes[n2])))
            }
          }
        }
      }
    }
    if (length(extendedConnections)>0) {
      adj_matAfter[startNode, names(extendedConnections)] = extendedConnections
      adj_matAfter[names(extendedConnections), startNode] = extendedConnections
    }
    adj_mat = adj_matAfter
  }
  startNodeNeighbors = G[names(which(abs(adj_mat[,startNode]) > 0))]
  startNodeUnvisitedNeighbors = startNodeNeighbors[!(names(startNodeNeighbors) %in% names(visitedNodes))]

  if (length(startNodeUnvisitedNeighbors)>0 || length(extendedConnections)>0) {
    weighted_edges.sum = sum(abs(adj_mat[which(rownames(adj_mat) %in% names(startNodeUnvisitedNeighbors)), startNode]))
    if (verbose==TRUE) {
      print(sprintf("%sWeighted_edges.sum=%f", paste(rep("   ", length(visitedNodes)-1), collapse=""), weighted_edges.sum))
    }
    z=1;
    while (z<=length(startNodeUnvisitedNeighbors)) {
      #Let startNodeUnvisitedNeighbors[z] be referred to as N from here on out.
      #N receives equal part of p1 amongst all of startNode's children.
      inherited.probability = p1*abs(adj_mat[names(startNodeUnvisitedNeighbors[z]),startNode])/weighted_edges.sum
      G[[names(startNodeUnvisitedNeighbors[z])]] = G[[names(startNodeUnvisitedNeighbors[z])]] + inherited.probability;

      if (verbose==TRUE) {
        print(sprintf("%sadded diffused  probability %f to child #%d: %s", paste(rep("   ", length(visitedNodes)-1), collapse=""), inherited.probability, z, names(startNodeUnvisitedNeighbors[z])));
      }

      # Determine if N has neighbors
      nNeighbors = G[names(which(abs(adj_mat[,names(startNodeUnvisitedNeighbors[z])])>0))]
      if (length(nNeighbors) > 0 && inherited.probability/2>thresholdDiff && ((length(visitedNodes)+1)<length(G))) {
        G[[names(startNodeUnvisitedNeighbors[z])]] = G[[names(startNodeUnvisitedNeighbors[z])]] - inherited.probability/2;
        if (verbose==TRUE) {
          print(sprintf("%ssubtracted %f from startNode's neighbor #%d: %s and sent to its own neighbors.",
                        paste(rep("   ", length(visitedNodes)-1), collapse=""), inherited.probability/2, z, names(startNodeUnvisitedNeighbors[z])));
        }
        G = graph.diffuseP1(inherited.probability/2, names(startNodeUnvisitedNeighbors[z]), G, c(visitedNodes, startNodeUnvisitedNeighbors[z]), graphNumber, verbose=FALSE);
      }
      z=z+1;
    }
  } else {
    # startNode is stranded with its visited neighbors, or is a singleton. Diffuse p1 uniformly amongst all unvisited nodes
    if (verbose==TRUE) {
      print("startNode is stranded with its visited neighbors, or is a singleton. Diffuse p1 uniformly amongst all unvisited nodes")
    }
    t=1;
    while (t<=length(G)) {
      if (!(names(G[t]) %in% names(visitedNodes))) {
        G[[t]] = G[[t]] + p1/(length(G)-length(visitedNodes));
      }
      t=t+1;
    }
  }
  return (G);
}


