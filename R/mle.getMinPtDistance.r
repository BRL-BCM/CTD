#' Get minimum patient distances
#' 
#' Given a series of patient distance matrices, return the minimum
#' distance between all pairwise patient comparisons made.
#' @param allSimMatrices - A list of all similarity matrices, across all 
#'                         k for a given graph, or across many graphs.
#' @return minPtSim - Pairwise patient distances representing the minimum
#' patient distance observed across several distance matrices.
#' @export mle.getMinPtDistance
#' @examples
#' # Get patient distances for the first 2 patients in the Miller 2015 dataset.
#' data("Miller2015")
#' data_mx = Miller2015[-c(1,grep("x - ",rownames(Miller2015))),
#'                         grep("IEM", colnames(Miller2015))]
#' data_mx = apply(data_mx[,c(1,2)], c(1,2), as.numeric)
#' # Build a network, G
#' adj_mat = matrix(0, nrow=nrow(data_mx), ncol=nrow(data_mx))
#' rows = sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' cols = sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' for(i in rows){for(j in cols){adj_mat[i,j]=rnorm(1,0,1)}}
#' colnames(adj_mat) = rownames(data_mx)
#' rownames(adj_mat) = rownames(data_mx)
#' G = vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat) 
#' # Look at the top 5 metabolites for each patient. 
#' kmx=5
#' topMets_allpts = c()
#' for(pt in seq_len(ncol(data_mx))){
#'     topMets_allpts=c(topMets_allpts, 
#'                     rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                             decreasing=TRUE)[seq_len(kmx)]])
#' }
#' topMets_allpts = unique(topMets_allpts)
#' # Pre-compute node ranks for all metabolites in topMets_allpts for
#' # faster distance calculations.
#' ranks = list()
#' for(n in seq_len(length(topMets_allpts))){ 
#'     ind=which(names(G)==topMets_allpts[n])
#'     ranks[[n]]=singleNode.getNodeRanksN(ind,G,0.9,0.01,adj_mat,
#'                                         topMets_allpts,log2(length(G))) 
#' }
#' names(ranks) = topMets_allpts
#' # Also pre-compute patient bitstrings for faster distance calculations.
#' ptBSbyK = list()
#' for (pt in seq_len(ncol(data_mx))) {
#'     S=rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                 decreasing=TRUE)[seq_len(kmx)]]
#'     ptBSbyK[[pt]]=mle.getPtBSbyK(S, ranks)
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
#'         tmp = mle.getPtDist(ptBSbyK[[pt]],ptID,ptBSbyK[[pt2]],ptID2,data_mx,
#'                             ranks,p1=0.9,thresholdDiff=0.01,adj_mat)
#'         for (k in seq_len(kmx)) {
#'             res[[k]]$ncd[ptID, ptID2] = tmp$NCD[k]
#'             res[[k]]$ncd[ptID2, ptID] = tmp$NCD[k]
#'         }
#'     }
#' }
#' res_ncd = lapply(res, function(i) i$ncd)
#' minPtDist = mle.getMinPtDistance(res_ncd)
mle.getMinPtDistance = function(allSimMatrices) {
    minPtSim = matrix(1000, nrow=nrow(allSimMatrices[[1]]),
                        ncol=ncol(allSimMatrices[[1]]))
    for (ind in seq_len(length(allSimMatrices))) {
        for (n1 in seq_len(nrow(allSimMatrices[[ind]]))) {
            for (n2 in seq_len(ncol(allSimMatrices[[ind]]))) {
                minPtSim[n1, n2] = min(minPtSim[n1, n2],
                                        allSimMatrices[[ind]][n1, n2])
            }
        }
    }
    return(minPtSim)
}