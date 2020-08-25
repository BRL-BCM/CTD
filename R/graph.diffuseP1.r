#' Diffuse Probability P1 from a starting node.
#'
#' Recursively diffuse probability from a starting node based on the connectivity of the background knowledge graph, 
#' representing the likelihood that a variable will be most influenced by a perturbation in the starting node.
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
#' @param recursion_level - The current depth in the call stack caused by a recursive algorithm. Only relevant if output_dir
#'                          is specified.
#' @param coords - The x and y coordinates for each node in the network, to remain static between images. Only relevant if
#'                 output_dir is specified.
#' @return G - A list of returned probabilities after the diffusion of probability has truncated, with names of the list 
#'             being the node names in the background knowledge graph.
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
#' probs_afterCurrDraw=graph.diffuseP1(p1=1.0, startNode=names(G)[1], G=G[1], visitedNodes=names(G)[1], 
#'                                     thresholdDiff=0.01, adj_mat, TRUE)
#' # Make a movie of the diffusion of probability from startNode
#' .GlobalEnv$imgNum = 1
#' probs_afterCurrDraw=graph.diffuseP1(p1=1.0, startNode=names(G)[1], G=G[1], visitedNodes=names(G)[1], 
#'                                     thresholdDiff=0.01, adj_mat, TRUE, getwd(), 1, 1)
graph.diffuseP1=function (p1, startNode, G, visitedNodes, thresholdDiff, adj_mat, verbose=FALSE, 
                          output_dir="", recursion_level=1, coords=NULL) {
  if (verbose==TRUE) {
    print(sprintf("%sprob. to diffuse:%f startNode: %s, visitedNodes: %s", 
                  paste(rep("   ",recursion_level-1),collapse=""),
                  p1, startNode, toString(visitedNodes)))
  }
  if (output_dir!="") { graph.takeDiffusionSnapShot(adj_mat, G, output_dir, p1, startNode, 
                                                    visitedNodes, recursion_level, coords) }
  adj_matAfter=graph.connectToExt(adj_mat, startNode, visitedNodes)
  startNodeNeighbors=names(which(abs(adj_matAfter[, startNode])>0))
  startNodeUnvisitedNeighbors=startNodeNeighbors[!(startNodeNeighbors %in% visitedNodes)]
  if (length(startNodeUnvisitedNeighbors)>0) {
    weighted_edges.sum=sum(abs(adj_matAfter[which(rownames(adj_matAfter) %in% startNodeUnvisitedNeighbors), startNode]))
    z=1
    while (z<=length(startNodeUnvisitedNeighbors)) {
      inherited.probability=p1*abs(adj_matAfter[startNodeUnvisitedNeighbors[z], startNode])/weighted_edges.sum
      G[[startNodeUnvisitedNeighbors[z]]]=G[[startNodeUnvisitedNeighbors[z]]]+inherited.probability
      if (verbose==TRUE) {
        print(sprintf("%sadded diffused  probability %f to child #%d: %s", 
                      paste(rep("   ",recursion_level-1),collapse=""), 
                      inherited.probability, z, startNodeUnvisitedNeighbors[z]))
      }
      if (output_dir!="") { graph.takeDiffusionSnapShot(adj_mat, G, output_dir, p1, startNode, 
                                                        visitedNodes, recursion_level, coords) }
      nNeighbors=G[names(which(abs(adj_matAfter[, startNodeUnvisitedNeighbors[z]])> 0))]
      if (length(nNeighbors)>0 && inherited.probability/2>thresholdDiff && ((length(visitedNodes) + 1) < length(G))) {
        G[[startNodeUnvisitedNeighbors[z]]]=G[[startNodeUnvisitedNeighbors[z]]]-inherited.probability/2
        if (verbose==TRUE) {
          print(sprintf("%ssubtracted %f from startNode's neighbor #%d: %s and sent to its own neighbors.",
                        paste(rep("   ", recursion_level-1),collapse=""), inherited.probability/2,
                        z, startNodeUnvisitedNeighbors[z]))
        }
        if (output_dir!="") { imgNum = graph.takeDiffusionSnapShot(adj_mat, G, output_dir, p1, startNode, 
                                                                   visitedNodes, recursion_level, coords) }
        
        G=graph.diffuseP1(inherited.probability/2, startNodeUnvisitedNeighbors[z], G, 
                          c(visitedNodes, startNodeUnvisitedNeighbors[z]), thresholdDiff, adj_mat,
                            verbose=verbose, output_dir, recursion_level+1, coords)
      }
      z=z + 1
    }
    if (output_dir!="") { graph.takeDiffusionSnapShot(adj_mat, G, output_dir, p1, startNode, 
                                                      visitedNodes, recursion_level, coords) }
  } else {
    if (verbose==TRUE) {
      print(sprintf("startNode %s is stranded with its visited neighbors, or is a singleton. Diffuse p1 uniformly amongst all unvisited nodes.", startNode))
    }
    G[!(names(G) %in% visitedNodes)]=unlist(G[!(names(G) %in% visitedNodes)]) + p1/(length(G) - length(visitedNodes))
    if (output_dir!="") { graph.takeDiffusionSnapShot(adj_mat, G, output_dir, p1, startNode, 
                                                      visitedNodes, recursion_level, coords) }
  }
  return(G)
}
