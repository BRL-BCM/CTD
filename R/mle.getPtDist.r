#' CTDncd: A network-based distance metric.
#'
#' This function calculates the universal distance between patients, using a
#' mutual information metric, where self-information comes from the minimum
#' encoding length of each patient's encoded modular perturbations in the
#' network.
#' @param p1.optBS - The optimal bitstring associated with patient 1.
#' @param ptID - The identifier associated with patient 1's sample.
#' @param p2.optBS - The optimal bitstring associated with patient 2.
#' @param ptID2 - The identifier associated with patient 2's sample.
#' @param data_mx - The matrix that gives the perturbation strength
#'                  (z-scores) for all variables (columns) for each
#'                  patient (rows).
#' @param ranks - The list of node ranks, starting with each node in patient
#'                1&2's subsets of interest.
#' @param p1 - The probability that is preferentially distributed between
#'             network nodes by the probability diffusion algorithm based
#'             solely on network connectivity. The remaining probability
#'             (i.e., "p0") is uniformally distributed between network nodes,
#'             regardless of connectivity.
#' @param thresholdDiff - When the probability diffusion algorithm exchanges
#'                        this amount (thresholdDiff) or less between nodes,
#'                        the algorithm returns up the call stack.
#' @param adj_mat - The adjacency matrix that encodes the edge weights for the
#'                  network, G. 
#' @return patientDistances - a distance matrix, where row and columns are
#' patient identifiers.
#' @export mle.getPtDist
#' @usage mle.getPtDist(p1.optBS,ptID,p2.optBS,ptID2,data_mx,ranks,p1,
#'                         thresholdDiff,adj_mat)
#' @examples
#' # Get patient distances for the first 2 patients in the Miller 2015 dataset.
#' data("Miller2015")
#' data_mx = Miller2015[-c(1,grep("x - ", rownames(Miller2015))),
#'                         grep("IEM",colnames(Miller2015))]
#' data_mx = apply(data_mx[,c(1,2)], c(1,2), as.numeric)
#' # Build a network, G
#' adj_mat = matrix(0, nrow=nrow(data_mx), ncol=nrow(data_mx))
#' rows = sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' cols = sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' for (i in rows) {for (j in cols) {adj_mat[i,j]=rnorm(1,mean=0,sd=1)}}
#' colnames(adj_mat) = rownames(data_mx)
#' rownames(adj_mat) = rownames(data_mx)
#' G = vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat) 
#' # Look at the top 5 metabolites for each patient. 
#' kmx=5
#' topMets_allpts = c()
#' for (pt in seq_len(ncol(data_mx))) {
#'     topMets_allpts=c(topMets_allpts, 
#'                     rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                             decreasing=TRUE)[seq_len(kmx)]])}
#' topMets_allpts = unique(topMets_allpts)
#' # Pre-compute node ranks for all metabolites in topMets_allpts
#' # for faster distance calculations.
#' ranks = list()
#' for (n in seq_len(length(topMets_allpts))) { 
#'     ind = which(names(G)==topMets_allpts[n])
#'     ranks[[n]]=singleNode.getNodeRanksN(ind,G,0.9,0.01,adj_mat,
#'                                         topMets_allpts,log2(length(G))) 
#' }
#' names(ranks) = topMets_allpts
#' # Also pre-compute patient bitstrings for faster distance calculations.
#' ptBSbyK = list()
#' for (pt in seq_len(ncol(data_mx))) {
#'     S=rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                 decreasing=TRUE)[seq_len(kmx)]]
#'     ptBSbyK[[pt]] = mle.getPtBSbyK(S, ranks)
#' }
#' # Build your results ("res") list object to store patient distances at
#' # different size k's.
#' res = list()
#' t = list(ncd=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
#' rownames(t$ncd) = colnames(data_mx)
#' colnames(t$ncd) = colnames(data_mx)
#' for (i in seq_len(kmx)) { res[[i]] = t }
#' for (pt in seq_len(ncol(data_mx))) {
#'     print(pt)
#'     ptID = colnames(data_mx)[pt]
#'     for (pt2 in pt:ncol(data_mx)) {
#'         ptID2 = colnames(data_mx)[pt2]
#'         tmp=mle.getPtDist(ptBSbyK[[pt]],ptID,ptBSbyK[[pt2]],ptID2,
#'                             data_mx,ranks,p1=0.9,thresholdDiff=0.01,adj_mat)
#'         for (k in seq_len(kmx)) {
#'             res[[k]]$ncd[ptID, ptID2] = tmp$NCD[k]
#'             res[[k]]$ncd[ptID2, ptID] = tmp$NCD[k]
#'         }
#'     }
#' }
mle.getPtDist=function(p1.optBS,ptID,p2.optBS,ptID2,data_mx,ranks,p1,
                        thresholdDiff,adj_mat) {
    if (length(p1.optBS)!=length(p2.optBS)){
        return("Error: Pt1 Subset different size from Pt2.")}
    if (is.null(ranks)) {return("Error: Must specify ranks")}
    G = vector(mode="list", length=nrow(data_mx))
    names(G) = rownames(data_mx)
    res.p1 = mle.getEncodingLength(p1.optBS, NULL, ptID, G)
    res.p2 = mle.getEncodingLength(p2.optBS, NULL, ptID, G)
    dirSim = stat.getDirSim(ptID, ptID2, length(p1.optBS), data_mx)
    #Get optimal encoding of patient1's union patient2's subsets
    p1.e=p2.e=p12.e=c()
    bits_fixed=log2(choose(length(G),1)) #fixed length encoding length
    for (k in seq_len(length(p1.optBS))) {
        p1.e[k]=res.p1[which.max(res.p1[seq_len(k),"d.score"]),"IS.alt"]+
            bits_fixed*(k-which.max(res.p1[seq_len(k),"d.score"]))
        p2.e[k]=res.p2[which.max(res.p2[seq_len(k),"d.score"]),"IS.alt"]+
            bits_fixed*(k-which.max(res.p2[seq_len(k),"d.score"]))
        # p12
        p1.sig.nodes_cpy=names(sort(abs(data_mx[,ptID]),
                                    decreasing=TRUE)[seq_len(length(p1.optBS))])
        p2.sig.nodes_cpy=names(sort(abs(data_mx[,ptID2]),
                                    decreasing=TRUE)[seq_len(length(p2.optBS))])
        p1.sig.nodes_k = names(which(p1.optBS[[k]]==1))
        p2.sig.nodes_k = names(which(p2.optBS[[k]]==1))
        while (length(p1.sig.nodes_k)<k) {
            p1.sig.nodes_k = unique(c(p1.sig.nodes_k, p1.sig.nodes_cpy[1]))
            p1.sig.nodes_cpy = p1.sig.nodes_cpy[-1]
        }
        while (length(p2.sig.nodes_k)<k) {
            p2.sig.nodes_k = unique(c(p2.sig.nodes_k, p2.sig.nodes_cpy[1]))
            p2.sig.nodes_cpy = p2.sig.nodes_cpy[-1]
        }
        p12.sig.nodes_k=vapply(unique(c(p1.sig.nodes_k, p2.sig.nodes_k)),
                                trimws, character(1))
        p12.optBS = mle.getPtBSbyK(p12.sig.nodes_k, ranks)
        res = mle.getEncodingLength(p12.optBS, NULL, NULL, G)
        p12.e[k] = res[which.max(res[,"d.score"]), "IS.alt"] +
            log2(choose(length(G), 1))*(length(p12.sig.nodes_k)-
                                        which.max(res[,"d.score"]))
    }
    # Normalized Compression Distance, Jaccard Distance (w/ Directionality)
    ncd=(p12.e-apply(cbind(p1.e,p2.e),1,min))/apply(cbind(p1.e,p2.e),1,max)
    ncd[which(ncd<0)]=0
    return(list(p1.e=p1.e, p2.e=p2.e, p12.e=p12.e, dirSim=dirSim, NCD=ncd))
}