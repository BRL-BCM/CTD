#' Capture the current location of a network walker
#' 
#' A network walker steps towards the node that inherited the highest
#' probability from the last node that it stepped into.
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
#' @param visitedNodes - A character vector of node names, storing the history
#'                       of previous draws in the node ranking.
#' @param S - A character vector of node names in the subset you want the
#'            network walker to find.
#' @param coords - The x and y coordinates for each node in the network, to 
#'                 remain static between images.
#' @param imgNum - The image number for this snapshot. If images are being 
#'                 generated in a sequence, this serves as an iterator for file
#'                 naming.
#' @param useLabels - If TRUE, node names will display next to their respective 
#'                    nodes in the network. If FALSE, node names will not
#'                    display.
#' @return 0
#' @export graph.netWalkSnapShot
#' @usage graph.netWalkSnapShot(adj_mat,G,output_dir,p1,visitedNodes,S,
#'                                 coords,imgNum=1,useLabels=TRUE)
#' @import igraph
#' @importFrom grDevices dev.off png
#' @importFrom graphics legend title
#' @examples
#' # 7 node example graph illustrating diffusion of probability based on network
#' # connectivity
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
#' G[seq_len(length(G))] = 0
#' names(G) = c("A", "B", "C", "D", "E", "F", "G")
#' S = c("A", "C")
#' coords = layout.fruchterman.reingold(ig)
#' # Uncomment to run
#' #graph.netWalkSnapShot(adj_mat,G,output_dir=getwd(),p1=1.0,
#' #                        "A",S,coords,1,TRUE)
graph.netWalkSnapShot = function(adj_mat, G, output_dir, p1, visitedNodes, S,
                                    coords, imgNum=1, useLabels=TRUE) {
    ig=graph.adjacency(adj_mat,mode="undirected",weighted=TRUE)
    V(ig)$color=rep("white", length(G))
    V(ig)$color[which(V(ig)$name %in% S)]="green"
    if(useLabels==TRUE){
        V(ig)$label = V(ig)$name
    }else{V(ig)$label=rep("", length(V(ig)$name))}
    png(sprintf("%s/netWalkMovie%d_%d.png", output_dir, 
                which(S==visitedNodes[1]), imgNum), 500, 500)
    plot.igraph(ig, layout=coords, vertex.color=V(ig)$color,
                vertex.label=V(ig)$label, edge.width=5*abs(E(ig)$weight),
                edge.arrow.size=0.01, vertex.size=5+round(50*unlist(G), 0), 
                mark.col="dark red", mark.groups = visitedNodes)
    title(sprintf("{%s}", paste(as.numeric(visitedNodes %in% S), collapse="")), 
            cex.main=2)
    legend("topright", legend=c("Jackpot Nodes", "Drawn Nodes"), 
            fill=c("green", "dark red"))
    dev.off()
    return(0)
}