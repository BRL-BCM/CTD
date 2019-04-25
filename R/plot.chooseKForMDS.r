#' Metabolite set enrichment analysis (MSEA) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' @param allSimMatrices - A list of all similarity matrices, across all k for a given graph, or across many graphs.
#' @param diagnoses - A character vector of diagnostic labels associated with the colnames of each matrix in allSimMatrices.
#' @param diseaseGraph - A diagnostic label found in diagnoses that the background knowledge graph was based on.
#' @export plot.chooseKForMDS.r
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
#' res_tmp = lapply(res, function(i) i$ncd)
#  aadc_k = plot.chooseKForMDS(res_tmp)
plot.chooseKForMDS = function(allSimMatrices, diagnoses, diseaseGraph) {
  names(diagnoses) = rownames(allSimMatrices[[1]])
  dist_cntds = vector("numeric", length=length(allSimMatrices))
  tpr = vector("numeric", length=length(allSimMatrices))
  tnr = vector("numeric", length=length(allSimMatrices))
  acc = vector("numeric", length=length(allSimMatrices))
  for (k in 1:length(allSimMatrices)) {
    patientSim = allSimMatrices[[k]]
    # Kmeans
    dis.centroid = apply(patientSim[which(diagnoses==diseaseGraph),], 2, mean)
    negCntl.centroid = apply(patientSim[which(diagnoses!=diseaseGraph),], 2, mean)
    kmns = kmeans(patientSim, centers = rbind(dis.centroid, negCntl.centroid))
    dcntd = diagnoses[names(kmns$cluster[which(kmns$cluster==1)])]
    ncntd = diagnoses[names(kmns$cluster[which(kmns$cluster==2)])]
    dist_cntds[k] = dist(rbind(dis.centroid, negCntl.centroid))
    tp = length(which(dcntd==diseaseGraph))
    fp = length(which(dcntd!=diseaseGraph))
    tn = length(which(ncntd!=diseaseGraph))
    fn = length(which(ncntd==diseaseGraph))
    tpr[k] = tp/(tp+fn)
    tnr[k] = tn/(tn+fp)
    acc[k] = (tp+tn)/(tp+tn+fn+fp)
    print(sprintf("K=%d, dist btn centroids = %.2f, TPR=%.2f, TNR=%.2f", k, dist_cntds[k], tpr[k], tnr[k]))
  }

  ind = which.max(dist_cntds+tpr+tnr)
  print(sprintf("K=%d selected, overall accuracy = %.2f.", ind, acc[ind]))


  return(ind)
}
