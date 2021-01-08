#' Connect a node to its unvisited "extended" neighbors
#' 
#' @param adj_mat - The adjacency matrix that encodes the edge weights 
#'                  for the network. 
#' @param startNode - The node most recently visited by the network walker,
#'                    from which p1 gets dispersed.
#' @param visitedNodes - The history of previous draws in the node ranking
#'                       sequence.
#' @return adj_matAfter - The adjacency matrix where the startNode is now
#' connected to its unvisited "extended" neighbors. An extended neighbor is
#' the neighbor of a neighbor.
#' @export graph.connectToExt
#' @examples
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
#' ig = graph.adjacency(as.matrix(adj_mat), mode="undirected",weighted=TRUE)
#' G=vector(mode="list", length=7)
#' G[seq_len(length(G))] = 0
#' names(G) = c("A", "B", "C", "D", "E", "F", "G")
#' startNode = "A"
#' visitedNodes = c("B", "C")
#' coords = layout.fruchterman.reingold(ig)
#' V(ig)$x = coords[,1]
#' V(ig)$y = coords[,2]
#' adj_matAfter = graph.connectToExt(adj_mat, startNode, visitedNodes)
graph.connectToExt=function(adj_mat, startNode, visitedNodes) {
    startNodeNbors=names(which(abs(adj_mat[, startNode])> 0))
    startNodeUnvisitedNbors=startNodeNbors[!(startNodeNbors%in%visitedNodes)]
    vN=visitedNodes[which(visitedNodes != startNode)]
    extConnections=NULL
    if (length(vN)>0 && length(startNodeUnvisitedNbors)==0) {
        adj_matAfter=adj_mat[-which(rownames(adj_mat) %in% vN),
                                -which(colnames(adj_mat) %in% vN)]
        connections=adj_mat[startNode, ]
        connYes=connections[which(abs(connections)>0)]
        connNo=connections[intersect(which(connections==0),
                                        which(!(names(connections)%in%
                                                    c(startNode, vN))))]
        if (length(connNo)>0) {
            for (n1 in seq_len(length(connNo))) {
                if (length(connYes)>0) {
                    for (n2 in seq_len(length(connYes))) {
                        if (abs(adj_mat[names(connYes[n2]),
                                        names(connNo[n1])])>0) {
                            connNo[n1]=adj_mat[names(connYes[n2]),
                                                names(connNo[n1])]
                            extConnections=c(extConnections,connNo[n1])
                        }
                    }
                }
            }
        }
        if (length(extConnections)>0) {
            adj_matAfter[startNode, names(extConnections)]=extConnections
            adj_matAfter[names(extConnections), startNode]=extConnections
        }
    } else {adj_matAfter = adj_mat}
    return(adj_matAfter)
}
