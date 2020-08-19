#' Diffuse Probability P1 from a starting node.
#'
#' Recursively diffuse probability from a starting node based on the connectivity of the background knowledge graph, representing the likelihood that a variable will be
#'         most influenced by a perturbation in the starting node.
#' @param p1 - The probability being dispersed from the starting node, startNode, which is preferentially distributed 
#'             between network nodes by the probability diffusion algorithm based solely on network connectivity.
#' @param startNode - The node most recently visited by the network walker, from which p1 gets dispersed.
#' @param G - A list of probabilities, with names of the list being the node names in the background knowledge graph.
#' @param visitedNodes - The history of previous draws in the node ranking sequence.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges this amount (thresholdDiff)
#'                        or less between nodes, the algorithm returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for the network, G. 
#' @param verbose - If debugging or tracking a diffusion event, verbose=TRUE will activate print statements. Default is FALSE.
#' @param output_dir - If specified, a image sequence will generate in the output directory specified.
#' @return G - A list of returned probabilities after the diffusion of probability has truncated, with names of the list being the node names in the background knowledge graph.
#' @export graph.diffuseP1
#' @keywords probability diffusion
#' @keywords network walker
#' @examples
#' # Read in any network via its adjacency matrix
#' adj_mat=matrix(1, nrow=100, ncol=100)
#' for (i in seq_len(100)) {for (j in seq_len(100)){adj_mat[i, j]=rnorm(1, mean=0, sd=1)} }
#' colnames(adj_mat)=sprintf("Metabolite%d", seq_len(100))
#' rownames(adj_mat)=colnames(adj_mat)
#' G=vector(mode="list", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' G=lapply(G, function(i) i[[1]]=0)
#' probs_afterCurrDraw=graph.diffuseP1(p1=1.0, startNode=names(G)[1], G=G[1], visitedNodes=names(G)[1], thresholdDiff=0.01, adj_mat, TRUE)
#' # Make a movie of the diffusion of probability from startNode
#' probs_afterCurrDraw=graph.diffuseP1(p1=1.0, startNode=names(G)[1], G=G[1], visitedNodes=names(G)[1], thresholdDiff=0.01, adj_mat, TRUE, getwd())
graph.diffuseP1=function (p1, startNode, G, visitedNodes, thresholdDiff, adj_mat, verbose=FALSE, output_dir="") {
  if (verbose==TRUE) {
    print(sprintf("%sprob. to diffuse:%f startNode: %s, visitedNodes: %s",
                  paste(rep("   ", length(visitedNodes) - 1), collapse=""),
                  p1, startNode, toString(visitedNodes)))
  }
  if (output_dir!="") { takeSnapShot(ig, output_dir, p1, startNode, recursion_level) }
  startNodeNeighbors=names(which(abs(adj_mat[, startNode])> 0))
  startNodeUnvisitedNeighbors=startNodeNeighbors[!(startNodeNeighbors %in% visitedNodes)]
  vN=visitedNodes[which(visitedNodes != startNode)]
  extendedConnections=NULL
  if (length(vN)>0 && length(startNodeUnvisitedNeighbors)==0) {
    adj_matAfter=adj_mat[-which(rownames(adj_mat) %in% vN), -which(colnames(adj_mat) %in% vN)]
    connections=adj_mat[startNode, ]
    connectionsYes=connections[which(abs(connections)>0)]
    connectionsNo=connections[intersect(which(connections==0), which(!(names(connections) %in% c(startNode, vN))))]
    if (length(connectionsNo)>0) {
      for (n1 in seq_len(length(connectionsNo))) {
        if (length(connectionsYes)>0) {
          for (n2 in seq_len(length(connectionsYes))) {
            if (abs(adj_mat[names(connectionsYes[n2]), names(connectionsNo[n1])])>0) {
              connectionsNo[n1]=adj_mat[names(connectionsYes[n2]), names(connectionsNo[n1])]
              extendedConnections=c(extendedConnections, connectionsNo[n1])
            }
          }
        }
      }
    }
    if (length(extendedConnections)>0) {
      adj_matAfter[startNode, names(extendedConnections)]=extendedConnections
      adj_matAfter[names(extendedConnections), startNode]=extendedConnections
    }
    adj_mat=adj_matAfter
  }
  startNodeNeighbors=names(which(abs(adj_mat[, startNode])>0))
  startNodeUnvisitedNeighbors=startNodeNeighbors[!(startNodeNeighbors %in% visitedNodes)]
  if (length(startNodeUnvisitedNeighbors)>0 || length(extendedConnections)>0) {
    weighted_edges.sum=sum(abs(adj_mat[which(rownames(adj_mat) %in% startNodeUnvisitedNeighbors), startNode]))
    if (verbose==TRUE) {print(sprintf("%sWeighted_edges.sum=%f", paste(rep("   ", length(visitedNodes) - 1), collapse=""), weighted_edges.sum))}
    z=1
    while (z<=length(startNodeUnvisitedNeighbors)) {
      inherited.probability=p1*abs(adj_mat[startNodeUnvisitedNeighbors[z], startNode])/weighted_edges.sum
      G[[startNodeUnvisitedNeighbors[z]]]=G[[startNodeUnvisitedNeighbors[z]]]+inherited.probability
      if (verbose==TRUE) {
        print(sprintf("%sadded diffused  probability %f to child #%d: %s", paste(rep("   ", length(visitedNodes) - 1), collapse=""), inherited.probability, z, startNodeUnvisitedNeighbors[z]))
      }
      nNeighbors=G[names(which(abs(adj_mat[, startNodeUnvisitedNeighbors[z]])> 0))]
      if (length(nNeighbors)>0 && inherited.probability/2>thresholdDiff && ((length(visitedNodes) + 1) < length(G))) {
        G[[startNodeUnvisitedNeighbors[z]]]=G[[startNodeUnvisitedNeighbors[z]]]-inherited.probability/2
        if (verbose==TRUE) {
          print(sprintf("%ssubtracted %f from startNode's neighbor #%d: %s and sent to its own neighbors.",
                        paste(rep("   ", length(visitedNodes) - 1), collapse=""), inherited.probability/2,
                        z, startNodeUnvisitedNeighbors[z]))
        }
        if (movie==TRUE) { takeSnapShot(ig, output_dir, p1, startNode, recursion_level) }
        
        G=graph.diffuseP1(inherited.probability/2, startNodeUnvisitedNeighbors[z], G, 
                          c(visitedNodes, startNodeUnvisitedNeighbors[z]), thresholdDiff, adj_mat,
                            verbose=FALSE)
      }
      z=z + 1
    }
    if (output_dir!="") { takeSnapShot(ig, output_dir, p1, startNode, recursion_level) }
  } else {
    if (verbose==TRUE) {print("startNode is stranded with its visited neighbors, or is a singleton. Diffuse p1 uniformly amongst all unvisited nodes")}
    G[!(names(G) %in% visitedNodes)]=unlist(G[!(names(G) %in% visitedNodes)]) + p1/(length(G) - length(visitedNodes))
  }
  return(G)
}
