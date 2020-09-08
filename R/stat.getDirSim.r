#' DirSim: The Jaccard distance with directionality incorporated.
#' 
#' @param ptID - The identifier associated with patient 1's sample.
#' @param ptID2 - The identifier associated with patient 2's sample.
#' @param kmx - The number of top perturbations to consider in distance
#'              calculation.
#' @param data_mx - The matrix that gives the perturbation strength
#'                  (z-scores) for all variables (columns) for each
#'                  patient (rows).
#' @return dirSim - a distance matrix, where row and columns are
#'                  patient identifiers.
#' @export stat.getDirSim
#' @examples
#' # Get patient distances for the first 2 patients in the Miller 2015 dataset.
#' data("Miller2015")
#' data_mx = Miller2015[-c(1,grep("x - ", rownames(Miller2015))),
#'                      grep("IEM",colnames(Miller2015))]
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
#' # Look at the top 15 metabolites for each patient. 
#' kmx=15
#' topMets_allpts = c()
#' for (pt in seq_len(ncol(data_mx))) {
#'   topMets_allpts=c(topMets_allpts, 
#'                    rownames(data_mx)[order(abs(data_mx[,pt]),
#'                                            decreasing=TRUE)[seq_len(kmx)]])}
#' topMets_allpts = unique(topMets_allpts)
#' # Pre-compute node ranks for all metabolites in topMets_allpts
#' # for faster distance calculations.
#' ranks = list()
#' for (n in seq_len(length(topMets_allpts))) { 
#'   ind = which(names(G)==topMets_allpts[n])
#'   ranks[[n]]=singleNode.getNodeRanksN(ind,G,0.9,0.01,adj_mat,
#'                                       topMets_allpts,log2(length(G))) 
#' }
#' names(ranks) = topMets_allpts
#' # Also pre-compute patient bitstrings for faster distance calculations.
#' ptBSbyK = list()
#' for (pt in seq_len(ncol(data_mx))) {
#'   S=rownames(data_mx)[order(abs(data_mx[,pt]),
#'                             decreasing=TRUE)[seq_len(kmx)]]
#'   ptBSbyK[[pt]] = singleNode.getPtBSbyK(S, ranks)
#' }
#' # Build your results ("res") list object to store patient distances at
#' # different size k's.
#' res = list()
#' t = list(dir=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
#' rownames(t$dir) = colnames(data_mx)
#' colnames(t$dir) = colnames(data_mx)
#' for (i in seq_len(kmx)) { res[[i]] = t }
#' for (pt in seq_len(ncol(data_mx))) {
#'   print(pt)
#'   ptID=colnames(data_mx)[pt]
#'   for (pt2 in pt:ncol(data_mx)) {
#'     ptID2=colnames(data_mx)[pt2]
#'     tmp=stat.getDirSim(ptID,ptID2,kmx,data_mx)
#'     for (k in seq_len(kmx)) {
#'       res[[k]]$dir[ptID, ptID2]=tmp[k]
#'       res[[k]]$dir[ptID2, ptID]=tmp[k]
#'     }
#'   }
#' }
stat.getDirSim = function(ptID, ptID2, kmx, data_mx) {
    # Get optimal bitstring for encoding of patient1's union patient2's subsets
    dirSim = vector("numeric", length=kmx)
    for (k in seq_len(kmx)) {
        p1.sig=rownames(data_mx)[order(abs(data_mx[,ptID]),
                                        decreasing=TRUE)][seq_len(k)]
        p2.sig=rownames(data_mx)[order(abs(data_mx[,ptID2]),
                                        decreasing=TRUE)][seq_len(k)]
        p1.dirs = data_mx[p1.sig, ptID]
        p1.dirs[which(!(p1.dirs>0))] = 0
        p1.dirs[which(p1.dirs>0)] = 1
        p2.dirs = data_mx[p2.sig, ptID2]
        p2.dirs[which(!(p2.dirs>0))] = 0
        p2.dirs[which(p2.dirs>0)] = 1
        p1.sig = sprintf("%s%d", p1.sig, p1.dirs)
        p2.sig = sprintf("%s%d", p2.sig, p2.dirs)
        p12.sig=union(p1.sig,p2.sig)
        dirSim[k]=1-(length(intersect(p1.sig,p2.sig))/length(p12.sig))
    }
    return(dirSim)
}