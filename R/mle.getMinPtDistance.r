#' Metabolite set enrichment analysis (MSEA) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' @param allSimMatrices - A list of all similarity matrices, across all k for a given graph, or across many graphs.
#' @export mle.getMinPtDistance.r
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Get patient distances
#' data_mx.pvals = apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
#' res = list()
#' t = list(ncd=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)),
#'          dir=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
#' rownames(t$ncd) = colnames(data_mx)
#' colnames(t$ncd) = colnames(data_mx)
#' rownames(t$dir) = colnames(data_mx)
#' colnames(t$dir) = colnames(data_mx)
#' for (i in 1:kmx) {
#'   res[[i]] = t
#' }
#' for (pt in 1:ncol(data_mx)) {
#'   print(pt)
#'   ptID = colnames(data_mx)[pt]
#'   for (pt2 in pt:ncol(data_mx)) {
#'     ptID2 = colnames(data_mx)[pt2]
#'     tmp = mle.getPtSim(ptBSbyK[[ptID]], ptID, ptBSbyK[[ptID2]], ptID2, data_mx, perms)
#'     for (k in 1:kmx) {
#'       res[[k]]$ncd[ptID, ptID2] = tmp$NCD[k]
#'       res[[k]]$dir[ptID, ptID2] = tmp$dirSim[k]
#'       res[[k]]$ncd[ptID2, ptID] = tmp$NCD[k]
#'       res[[k]]$dir[ptID2, ptID] = tmp$dirSim[k]
#'     }
#'   }
#' }
#' patientSimilarity = mle.getMinPtDistance(res)
mle.getMinPtDistance = function(allSimMatrices) {
  ptSim = matrix(1000, nrow=nrow(allSimMatrices[[1]]), ncol=ncol(allSimMatrices[[1]]))
  for (ind in 1:length(allSimMatrices)) {
    for (n1 in 1:nrow(allSimMatrices[[ind]])) {
      for (n2 in 1:ncol(allSimMatrices[[ind]])) {
        ptSim[n1, n2] = min(ptSim[n1, n2], allSimMatrices[[ind]][n1, n2])
      }
    }
  }
  
  return(ptSim)
}
