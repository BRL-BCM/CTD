#' Generate single-node node rankings ("fixed" walk) 
#'
#' This function calculates the node rankings starting from a given perturbed 
#' variable in a subset of variables in the network.
#' @param n - The index (out of a vector of node names) of the node ranking
#'            you want to calculate.
#' @param G - A list of probabilities with list names being the node names
#'            of the network.
#' @param S - A character vector of node names in the subset you want the 
#'            network walker to find.
#' @param num.misses - The number of "misses" the network walker will tolerate
#'                     before switching to fixed length codes for remaining
#'                     nodes to be found.
#' @param p1 - The probability that is preferentially distributed between
#'             network nodes by the probability diffusion algorithm based
#'             solely on network connectivity. The remaining probability
#'             (i.e., "p0") is uniformally distributed between network nodes,
#'             regardless of connectivity.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges 
#'                        this amount or less between nodes, the algorithm 
#'                        returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for the 
#'                  network, G. 
#' @param verbose - If TRUE, print statements will execute as progress is made.
#'                  Default is FALSE.
#' @param out_dir - If specified, a image sequence will generate in the 
#'                  output directory specified.
#' @param useLabels - If TRUE, node names will display next to their respective
#'                    nodes in the network. If FALSE, node names will not
#'                    display. Only relevant if out_dir is specified. 
#' @param coords - The x and y coordinates for each node in the network, to 
#'                 remain static between images.
#' @return curr_ns - A character vector of node names in the order they
#' were drawn by the probability diffusion algorithm.
#' @keywords probability diffusion
#' @keywords network walker
#' @export singleNode.getNodeRanksN
#' @usage singleNode.getNodeRanksN(n,G,p1,thresholdDiff,adj_mat,
#'                                     S=NULL,num.misses=NULL,verbose=FALSE,
#'                                     out_dir="",useLabels=FALSE,coords=NULL)
#' @examples
#' # Build an adjacency matrix for network G
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
#' G=vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' # Get node rankings for the first metabolite in network G. 
#' ranks=singleNode.getNodeRanksN(1,G,p1=0.9,thresholdDiff=0.01,adj_mat)
singleNode.getNodeRanksN = function(n,G,p1,thresholdDiff,adj_mat,S=NULL,
                                    num.misses=NULL,verbose=FALSE,out_dir="",
                                    useLabels=FALSE,coords=NULL) {
    p0=1-p1
    if (is.null(S) && (!is.null(num.misses) || out_dir!="")) {
        print("You must also supply S if out_dir or num.misses is supplied")
        return(0)}
    if(verbose){print(sprintf("Node ranking %d of %d.",n,length(G)))}
    curr_ns = NULL # current node set
    stopIterating=FALSE
    startNode = names(G)[n]
    currGph = G
    numMisses = 0
    curr_ns = c(curr_ns, startNode)
    if(out_dir!=""){graph.netWalkSnapShot(adj_mat,G,out_dir,p1,curr_ns,S, 
                                            coords,length(curr_ns),useLabels)}
    while (stopIterating==FALSE) {
        currGph[seq_len(length(currGph))]=0 # clear probabilities
        baseP=p0/(length(currGph)-length(curr_ns))
        #set unvisited nodes to baseP
        currGph[!(names(currGph) %in% curr_ns)]=baseP
        currGph=graph.diffuseP1(p1,startNode,currGph,curr_ns,thresholdDiff,
                                adj_mat,verbose=FALSE)
        # Sanity check. p1_event should add up to roughly p1
        p1_event = sum(unlist(currGph[!(names(currGph) %in% curr_ns)]))
        if (abs(p1_event-1)>thresholdDiff) {
            extra.prob.to.diffuse=1-p1_event
            currGph[names(curr_ns)]=0
            ind=!(names(currGph)%in%names(curr_ns))
            currGph[ind]=unlist(currGph[ind])+extra.prob.to.diffuse/sum(ind)}
        # Set startNode to the node with the max probability in the new currGph
        maxProb=names(which.max(currGph))
        if(out_dir!=""){graph.netWalkSnapShot(adj_mat,G,out_dir,p1,curr_ns,S, 
                                                coords,length(curr_ns),
                                                useLabels)}
        # Break ties: When there are ties, choose the first of the winners.
        startNode = names(currGph[maxProb[1]])
        if (!is.null(S)) { # draw until all members of S are found
            if(startNode %in% S){numMisses=0}else{numMisses=numMisses+1}
            curr_ns = c(curr_ns, startNode)
            if (numMisses>num.misses||all(S %in% curr_ns)){stopIterating=TRUE}
        } else { # keep drawing until you've drawn all nodes in G
            curr_ns = c(curr_ns, startNode)
            if(length(curr_ns)>=(length(G))){stopIterating=TRUE}}
        if(out_dir!=""){graph.netWalkSnapShot(adj_mat,G,out_dir,p1,curr_ns,S, 
                                                coords,length(curr_ns),
                                                useLabels)}
    }
    return(curr_ns)
}