#' Generate multi-node node rankings ("adaptive" walk) 
#'
#' This function calculates the node rankings starting from a given node in a
#' subset of nodes in a given network, G.
#' @param S - A character vector of the node names for the subset of nodes you
#'            want to encode.
#' @param G - A list of probabilities with list names being the node names of
#'            the network.
#' @param p1 - The probability that is preferentially distributed between
#'             network nodes by the probability diffusion algorithm based
#'             solely on network connectivity. The remaining probability, 1-p1,
#'             is uniformally distributed between network nodes, regardless of
#'             connectivity.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges
#'                        this amount or less between nodes, the algorithm 
#'                        returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for the
#'                  network, G. 
#' @param num.misses - The number of "misses" the network walker will tolerate
#'                     before switching to fixed length codes for remaining
#'                     nodes to be found.
#' @param verbose - If TRUE, print statements will execute as progress is made.
#'                  Default is FALSE.
#' @param out_dir - If specified, a image sequence will generate in the
#'                     output directory specified.
#' @param useLabels - If TRUE, node names will display next to their respective
#'                    nodes in the network. If FALSE, node names will not
#'                    display. Only relevant if out_dir is specified. 
#' @param coords - The x and y coordinates for each node in the network, to
#'                 remain static between images.
#' @return ranks - A list of character vectors of node names in the order they
#' were drawn by the probability diffusion algorithm, from each starting node
#' in S.
#' @keywords probability diffusion algorithm
#' @keywords network walker algorithm
#' @export multiNode.getNodeRanks
#' @usage multiNode.getNodeRanks(S,G,p1,thresholdDiff,adj_mat,num.misses=NULL,
#'                                 verbose=FALSE,out_dir="",useLabels=FALSE,
#'                                 coords=NULL)
#' @examples
#' # Read in any network via its adjacency matrix
#' adj_mat=rbind(c(0,1,2,0,0,0,0,0,0), #A's neighbors
#'                 c(1,0,3,0,0,0,0,0,0), #B's neighbors
#'                 c(2,3,0,0,1,0,0,0,0), #C's neighbors
#'                 c(0,0,0,0,0,0,1,1,0), #D's neighbors
#'                 c(0,0,1,0,0,1,0,0,0), #E's neighbors
#'                 c(0,0,0,0,1,0,0,0,0), #F's neighbors
#'                 c(0,0,0,1,0,0,0,1,0), #G's neighbors
#'                 c(0,0,0,1,0,0,1,0,0), #H's neighbors
#'                 c(0,0,0,0,0,0,0,0,0) #I's neighbors
#'                 )
#' rownames(adj_mat)=c("A","B","C","D","E","F","G","H","I")
#' colnames(adj_mat)=c("A","B","C","D","E","F","G","H","I")
#' G=vector(mode="list", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' S=names(G)[seq_len(3)]
#' ranks=multiNode.getNodeRanks(S, G, p1=0.9, thresholdDiff=0.01, adj_mat)
multiNode.getNodeRanks=function(S,G,p1,thresholdDiff,adj_mat,num.misses=NULL,
                                verbose=FALSE,out_dir="",useLabels=FALSE,
                                coords=NULL){
    if(is.null(num.misses)){num.misses=log2(length(G))}
    ranks=list()
    for (n in seq_len(length(S))) {
        if(verbose){print(sprintf("Node rankings %d of %d.", n, length(S)))}
        stopIterating=FALSE
        curr_ns=startNode=S[n] #current node set
        numMisses=0
        if(out_dir!=""){graph.netWalkSnapShot(adj_mat,G,out_dir,p1,curr_ns,S, 
                                                coords,length(curr_ns),
                                                useLabels)}
        while(stopIterating==FALSE) {
            hits = curr_ns[which(curr_ns %in% S)]#diffuse from all nodes in hits
            sumHits=as.vector(matrix(0, nrow=length(G), ncol=1))
            names(sumHits)=names(G)
            for (hit in seq_len(length(hits))) {
                G[seq_len(length(G))]=0 #Clear probabilities.
                baseP=(1-p1)/(length(G)-length(curr_ns))
                G[which(!(names(G) %in% curr_ns))]=baseP
                G=graph.diffuseP1(p1,hits[hit],G,curr_ns,thresholdDiff,adj_mat)
                #Sanity check: p1_event should be within 'thresholdDiff' of p1.
                p1_event=sum(unlist(G[!(names(G) %in% curr_ns)]))
                if (abs(p1_event-1)>thresholdDiff) {
                    G[names(curr_ns)]=0
                    ind=which(!(names(G) %in% names(curr_ns)))
                    G[ind]=unlist(G[ind]) + (1-p1_event)/length(ind)}
                sumHits=sumHits+unlist(G)
            }
            sumHits=sumHits/length(hits) #Get avg prob across diffusion events
            maxProb=names(which.max(sumHits[-which(names(sumHits)%in%curr_ns)]))
            if(out_dir!=""){graph.netWalkSnapShot(adj_mat,G,out_dir,p1,curr_ns,
                                                    S,coords,length(curr_ns),
                                                    useLabels)}
            startNode=names(G[maxProb[1]]) #Break ties: Choose first winner
            if (startNode %in% S) {
                numMisses=0
                hits=c(hits, startNode)} else {numMisses=numMisses+1}
            curr_ns=c(curr_ns, startNode)
            if(all(S %in% curr_ns) || numMisses>num.misses){stopIterating=TRUE}
            if(out_dir!=""){graph.netWalkSnapShot(adj_mat,G,out_dir,p1,curr_ns,
                                                    S,coords,length(curr_ns),
                                                    useLabels)}
        }
        ranks[[n]]=curr_ns
    }
    names(ranks)=S
    return(ranks)
}