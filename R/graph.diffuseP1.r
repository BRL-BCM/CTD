#' Diffuse Probability P1 from a starting node
#'
#' Recursively diffuse probability from a starting node based on the
#' connectivity of the network, representing the likelihood that a 
#' variable is most influenced by a perturbation in the starting node.
#' @param p1 - The probability being dispersed from the starting node, 
#'             sn, which is preferentially distributed between 
#'             network nodes by the probability diffusion algorithm 
#'             based solely on network connectivity.
#' @param sn - "Start node", or the node most recently visited by the
#'             network walker, from which p1 gets dispersed.
#' @param G - A list of probabilities, with names of the list being the 
#'            node names in the network.
#' @param vNodes - "Visited nodes", or the history of previous draws
#'                 in the node ranking sequence.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges
#'                        this amount (thresholdDiff) or less between nodes, 
#'                        the algorithm returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for
#'                  the network, G. 
#' @param verbose - If debugging or tracking a diffusion event, verbose=TRUE
#'                  will activate print statements. Default is FALSE.
#' @param out_dir - If specified, a image sequence will generate in the
#'                  output directory specified.
#' @param r_level - "Recursion level", or the current depth in the call stack 
#'                  caused by a recursive algorithm. Only relevant if out_dir
#'                  is specified.
#' @param coords - The x and y coordinates for each node in the network, to
#'                 remain static between images. Only relevant if out_dir
#'                 is specified.
#' @return G - A list of returned probabilities after the diffusion of
#' probability has truncated, with names of the list being the node names
#' in the network.
#' @export graph.diffuseP1
#' @keywords probability diffusion
#' @keywords network walker
#' @usage graph.diffuseP1(p1,sn,G,vNodes,thresholdDiff,adj_mat,verbose=FALSE,
#'                         out_dir="",r_level=1,coords=NULL)
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
#' G=lapply(G, function(i) i[[1]]=0)
#' probs_afterCurrDraw=graph.diffuseP1(p1=1.0, sn=names(G)[1], G=G,
#'                                     vNodes=names(G)[1], 
#'                                     thresholdDiff=0.01, adj_mat, TRUE)
graph.diffuseP1=function(p1,sn,G,vNodes,thresholdDiff,adj_mat,verbose=FALSE,
                            out_dir="",r_level=1,coords=NULL){
    nTabs=paste(rep("   ",r_level-1),collapse="") # for verbose stmts
    if(verbose==TRUE) {
        print(sprintf("%sprob. to diffuse:%f sn: %s, visitedNodes: %s",
                        nTabs,p1,sn,toString(vNodes)))}
    if(out_dir!=""){graph.diffusionSnapShot(adj_mat,G,out_dir,p1,sn,
                                            vNodes,r_level,coords)}
    adj_mat2=graph.connectToExt(adj_mat,sn,vNodes)
    snNbors=names(which(abs(adj_mat2[,sn])>0)) # sn's neighbors
    snUvNbors=snNbors[!(snNbors%in%vNodes)] # sn's unvisited neighbors
    if(length(snUvNbors)>0) {
        wEgdes.sum=sum(abs(adj_mat2[which(rownames(adj_mat2)%in%snUvNbors),sn]))
        z=1
        while(z<=length(snUvNbors)) {
            i.prob=p1*abs(adj_mat2[snUvNbors[z],sn])/wEgdes.sum # inherited prob
            G[[snUvNbors[z]]]=G[[snUvNbors[z]]]+i.prob
            if(verbose==TRUE) {
            print(sprintf("%schild#%d %s got %f",nTabs,z,snUvNbors[z],i.prob))}
            if(out_dir!=""){graph.diffusionSnapShot(adj_mat,G,out_dir,p1,sn,
                                                    vNodes,r_level,coords)}
            nNbors=G[names(which(abs(adj_mat2[,snUvNbors[z]])>0))]
            if(length(nNbors)>0 && i.prob/2>thresholdDiff &&
                ((length(vNodes)+1)<length(G))) {
                    G[[snUvNbors[z]]]=G[[snUvNbors[z]]]-i.prob/2
                    if(verbose==TRUE){
                        print(sprintf("%stook %f from child#%d:%s to send",
                                        nTabs,i.prob/2,z,snUvNbors[z]))}
                    if(out_dir!=""){graph.diffusionSnapShot(adj_mat,G,out_dir,
                                                            p1,sn,vNodes,
                                                            r_level,coords)}
                    G=graph.diffuseP1(i.prob/2,snUvNbors[z],G,
                                        c(vNodes,snUvNbors[z]),thresholdDiff,
                                        adj_mat,verbose=verbose,out_dir,
                                        r_level+1,coords)}
            z=z + 1
        }
        if(out_dir!=""){graph.diffusionSnapShot(adj_mat,G,out_dir,p1,sn,
                                                vNodes,r_level,coords)}
    } else {
        if(verbose==TRUE) {
            print(sprintf("%s is singleton or stranded by visited n.bors",sn))
            print("Diffuse p1 uniformly amongst all unvisited nodes.")}
        G[!(names(G)%in%vNodes)]=unlist(G[!(names(G)%in%vNodes)])+
                                        p1/(length(G)-length(vNodes))
        if(out_dir!=""){graph.diffusionSnapShot(adj_mat,G,out_dir,p1,sn,
                                                vNodes,r_level,coords)}
    }
    return(G)
}
