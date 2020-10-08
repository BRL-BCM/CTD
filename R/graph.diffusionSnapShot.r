#' Capture the current state of probability diffusion
#'
#' Recursively diffuse probability from a starting node based on the
#' connectivity in a network, G, where the probability represents the
#' likelihood that a variable will be influenced by a perturbation
#' in the starting node.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for
#'                  the network, G. 
#' @param G - A list of probabilities, with names of the list being the node
#'            names in the network.
#' @param output_dir - The local directory at which you want still PNG images
#'                     to be saved.
#' @param p1 - The probability being dispersed from the starting node,
#'             startNode, which is preferentially distributed between network
#'             nodes by the probability diffusion algorithm based solely on
#'             network connectivity.
#' @param startNode - The first variable drawn in the node ranking, from
#'                    which p1 gets dispersed.
#' @param visitedNodes - A character vector of node names, storing the
#'                       history of previous draws in the node ranking.
#' @param recursion_level - The current depth in the call stack caused by
#'                          a recursive algorithm.
#' @param coords - The x and y coordinates for each node in the network, to
#'                 remain static between images.
#' @return 0
#' @export graph.diffusionSnapShot
#' @usage graph.diffusionSnapShot(adj_mat,G,output_dir,p1,startNode,
#'                                 visitedNodes,recursion_level,coords)
#' @import igraph
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend title
#' @examples
#' # 7 node example graph illustrating diffusion of probability based on
#' # network connectivity.
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
#' ig = graph.adjacency(as.matrix(adj_mat),mode="undirected",weighted=TRUE)
#' G=vector(mode="list", length=7)
#' G[seq_len(length(G))] = 0
#' names(G) = c("A", "B", "C", "D", "E", "F", "G")
#' coords = layout.fruchterman.reingold(ig)
#' V(ig)$x = coords[,1]
#' V(ig)$y = coords[,2]
#' # Uncomment to run
#' #graph.diffusionSnapShot(adj_mat,G,getwd(),1.0,"A","A",1,coords)
graph.diffusionSnapShot = function(adj_mat, G, output_dir, p1,
                                        startNode, visitedNodes,
                                        recursion_level=1, coords) {
    ig = graph.adjacency(adj_mat, mode="undirected", weighted = TRUE)
    G = G[which(names(G) %in% V(ig)$name)]
    V(ig)$color = rep("blue", length(G))
    V(ig)$color[which(V(ig)$name %in% visitedNodes)] = "red"
    vals = rep(0, length(V(ig)$name))
    names(vals) = V(ig)$name
    vals[which(names(vals) %in% names(G))] = G
    V(ig)$label = sprintf("%s:%.2f", V(ig)$name, vals)
    curr_time = unclass(as.POSIXlt(Sys.time()))
    curr_time = sprintf("%s-%s-%s-%s-%s-%.3f", curr_time$year, curr_time$mon, 
                        curr_time$mday, curr_time$hour, curr_time$min, 
                        as.numeric(curr_time$sec))
    png(sprintf("%s/diffusionP1Movie%s.png", output_dir, curr_time),500,500)
    plot.igraph(ig, layout=coords, vertex.color=V(ig)$color,
                vertex.label=V(ig)$label, vertex.label.dist = 3,
                edge.width=5*abs(E(ig)$weight), mark.col="black", 
                mark.border = "black", mark.groups = startNode)
    title(sprintf("Diffuse %.2f from %s at recursion level %d.",
                    p1, startNode, recursion_level), cex.main=1)
    legend("bottomright", legend=c("Visited", "Unvisited"),
            fill=c("red", "blue"))
    dev.off()
    return(0)
}


