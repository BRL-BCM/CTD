#' Generate patient-specific bitstrings
#'
#' This function calculates the bitstrings (1 is a hit; 0 is a miss)
#' associated with a network walker which tries to find all nodes in
#' a given subset, S, in a given network, G.
#' @param S - A character vector of node names describing the node subset
#'            to be encoded.
#' @param ranks - The list of node ranks calculated over all possible nodes,
#'                starting with each node in subset of interest.
#' @param num.misses - The number of misses tolerated by the network walker
#'                     before path truncation occurs.
#' @return pt.byK - a list of bitstrings, with the names of the list elements
#' the node names of the encoded nodes
#' @export mle.getPtBSbyK
#' @examples
#' # Get patient bitstrings for the first 2 patients in the Miller 2015 dataset.
#' data("Miller2015")
#' data_mx=Miller2015[-c(1,grep("x - ", rownames(Miller2015))),
#'                     grep("IEM", colnames(Miller2015))]
#' data_mx=apply(data_mx[,c(1,2)], c(1,2), as.numeric)
#' # Build an adjacency matrix for network G
#' adj_mat=matrix(0, nrow=nrow(data_mx), ncol=nrow(data_mx))
#' rows=sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' cols=sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' for(i in rows){for (j in cols){adj_mat[i, j]=rnorm(1,0,1)}}
#' colnames(adj_mat)=rownames(data_mx)
#' rownames(adj_mat)=rownames(data_mx)
#' G=vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' # Look at the top 5 metabolites for each patient. 
#' kmx=5
#' topMets_allpts=c()
#' for (pt in seq_len(ncol(data_mx))) { 
#'     topMets_allpts=c(topMets_allpts,
#'                     rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                             decreasing=TRUE)[seq_len(kmx)]])}
#' topMets_allpts=unique(topMets_allpts)
#' # Use a single-node or multi-node network walker.
#' # Here we use a single-node network walker.
#' ranks=list()
#' for (n in seq_len(length(topMets_allpts))) { 
#'     ind=which(names(G)==topMets_allpts[n])
#'     ranks[[n]]=singleNode.getNodeRanksN(ind,G,0.9,0.01,adj_mat,
#'                                         topMets_allpts,log2(length(G))) 
#' }
#' names(ranks)=topMets_allpts
#' ptBSbyK=list()
#' for (pt in seq_len(ncol(data_mx))) {
#'     S=rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                 decreasing=TRUE)[seq_len(kmx)]]
#'     ptBSbyK[[pt]]=mle.getPtBSbyK(S, ranks)
#' }
mle.getPtBSbyK=function(S, ranks, num.misses=NULL) {
    if (is.null(num.misses)) { num.misses=log2(length(ranks)) }
    ranks2=ranks[which(names(ranks) %in% S)]
    pt.bitString=list()
    for (p in seq_len(length(S))) {
        miss=0
        thresh=length(ranks2[[S[p]]])
        for (ii in seq_len(length(ranks2[[S[p]]]))) {
            ind_t=as.numeric(ranks2[[S[p]]][ii] %in% S)
            if (ind_t==0) {
                miss=miss + 1
                if (miss >= num.misses) {
                    thresh=ii
                    break;
                }
            } else { miss=0 }
            pt.bitString[[S[p]]][ii]=ind_t
        }
        pt.bitString[[S[p]]]=pt.bitString[[S[p]]][seq_len(thresh)]
        names(pt.bitString[[S[p]]])=ranks2[[S[p]]][seq_len(thresh)]
        ind=which(pt.bitString[[S[p]]] == 1)
        pt.bitString[[S[p]]]=pt.bitString[[S[p]]][seq_len(ind[length(ind)])]
    }
    # For each k, find the ranks that found at least k in S. Which node 
    # rankings, from different start nodes, found the first k soonest?
    pt.byK=list()
    for (k in seq_len(length(S))) {
        pt.byK_tmp=pt.bitString
        # Which found at least k nodes
        bestInd=vector("numeric", length(S))
        for (p in seq_len(length(S))) {
            bestInd[p]=sum(pt.bitString[[p]])
        }
        if (length(which(bestInd>=k))>0) {
            pt.byK_tmp=pt.byK_tmp[which(bestInd>=k)]
            # Which found the most nodes soonest
            bestInd=vector("numeric", length(pt.byK_tmp))
            for (p in seq_len(length(pt.byK_tmp))) {
                bestInd[p]=sum(which(pt.byK_tmp[[p]] == 1)[seq_len(k)])
            }
            maxInd=max(which(pt.byK_tmp[[which.min(bestInd)]]==1)[seq_len(k)])
            pt.byK[[k]]=pt.byK_tmp[[which.min(bestInd)]][seq_len(maxInd)]
        } else {pt.byK[[k]]=pt.byK_tmp[[which.max(bestInd)]]}
    }
    return(pt.byK)
}