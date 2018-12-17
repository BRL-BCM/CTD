#' Generate heatmap plot of patient similarity matrix.
#'
#' This function plots a heatmap of a patient similarity matrix.
#' @param simMat - The patient similarity matrix.
#' @param path - The filepath to a directory in which you want to store the .png file.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of simMat.
#' @export plot.hmSim
#' @examples
#' # Read in any network via its adjacency matrix
#' tmp = matrix(1, nrow=100, ncol=100)
#' for (i in 1:100) {
#'   for (j in 1:100) {
#'     tmp[i, j] = rnorm(1, mean=0, sd=1)
#'   }
#' }
#' colnames(tmp) = sprintf("MolPheno%d", 1:100)
#' ig = graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))  # Must have this declared as a GLOBAL variable!!!!!
#' # Set other tuning parameters
#' p0=0.1  # 10% of probability distributed uniformly
#' p1=0.9  # 90% of probability diffused based on edge weights in networks
#' thresholdDiff=0.01
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = V(ig)$name
#' # Get node permutations for graph
#' perms = list()
#' for (n in 1:length(G)) {
#'   print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
#'   perms[[n]] = mle.getPermN(n, G)
#' }
#' names(perms) = names(G)
#' # Decide what the largest subset size you will consider will be
#' kmx = 20
#' # Load your patient data (p features as rows x n observations as columns)
#' # data_mx = read.table("/your/own/data.txt", sep="\t", header=TRUE)
#' data(testData)
#' data_mx = t(testData)
#' rownames(data_mx) = tolower(rownames(data_mx))
#' # Get bitstrings associated with each patient's top kmx variable subsets
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'   ptID = colnames(data_mx)[pt]
#'   ptBSbyK[[ptID]] = mle.getPtBSbyK(data_mx, ptID, perms, kmx)
#' }
#' # Get patient distances
#' data_mx.pvals = apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
#' res = list()
#' t = list(ncd=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)),
#'          dir=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)),
#'          jac=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
#' rownames(t$ncd) = colnames(data_mx)
#' colnames(t$ncd) = colnames(data_mx)
#' rownames(t$dir) = colnames(data_mx)
#' colnames(t$dir) = colnames(data_mx)
#' rownames(t$jac) = colnames(data_mx)
#' colnames(t$jac) = colnames(data_mx)
#' for (i in 1:(kmx-1)) {
#'   res[[i]] = t
#' }
#' for (pt in 1:ncol(data_mx)) {
#'   print(pt)
#'   ptID = colnames(data_mx)[pt]
#'   for (pt2 in pt:ncol(data_mx)) {
#'     ptID2 = colnames(data_mx)[pt2]
#'     for (k in 1:(kmx-1)) {
#'       tmp = mle.getPatientSimilarity(ptBSbyK[[ptID]][k], ptID, ptBSbyK[[ptID2]][k], ptID2, data_mx, perms)
#'       res[[k]]$ncd[ptID, ptID2] = tmp$NCD
#'       res[[k]]$dir[ptID, ptID2] = tmp$dirSim
#'       res[[k]]$ncd[ptID2, ptID] = tmp$NCD
#'       res[[k]]$dir[ptID2, ptID] = tmp$dirSim
#'
#'       p1.sig.nodes = rownames(data_mx)[order(abs(data_mx[,ptID]), decreasing = TRUE)][1:k]
#'       p2.sig.nodes = rownames(data_mx)[order(abs(data_mx[,ptID2]), decreasing = TRUE)][1:k]
#'       p1.dirs = data_mx[p1.sig.nodes, ptID]
#'       p1.dirs[which(!(p1.dirs>0))] = 0
#'       p1.dirs[which(p1.dirs>0)] = 1
#'       p2.dirs = data_mx[p2.sig.nodes, ptID2]
#'       p2.dirs[which(!(p2.dirs>0))] = 0
#'       p2.dirs[which(p2.dirs>0)] = 1
#'       p1.sig.nodes = sprintf("%s%d", p1.sig.nodes, p1.dirs)
#'       p2.sig.nodes = sprintf("%s%d", p2.sig.nodes, p2.dirs)
#'       res[[k]]$jac[ptID, ptID2] = 1 - (length(intersect(p1.sig.nodes, p2.sig.nodes))/length(union(p1.sig.nodes, p2.sig.nodes)))
#'       res[[k]]$jac[ptID2, ptID] = 1 - (length(intersect(p1.sig.nodes, p2.sig.nodes))/length(union(p1.sig.nodes, p2.sig.nodes)))
#'     }
#'   }
#' }
#' # if you have diagnostic labels associated with the colnames(data_mx), send them using diagnoses parameter
#' diagnoses = colnames(data_mx)
#' diagnoses[1:25] = "diseased"
#' diagnoses[26:50] = "neg_control"
#' patientSim = 0.8*res[[k]]$ncd + 0.2*res[[k]]$dir
#' plot.hmSim(patientSim, path=getwd(), diagnoses)
plot.hmSim = function(simMat, path, diagnoses=NULL) {
  if (is.null(diagnoses)) {
    png(sprintf("%s/ptSimilarity.png", path), 3000, 1000)
    heatmap.2(x=simMat,
              dendrogram="both",
              Rowv=TRUE,
              Colv=TRUE,
              cexRow=0.75,cexCol=1, margins=c(12,12),
              trace="none", key=TRUE,
              col=bluered,
              notecol="black")
    dev.off()
  } else {
    d = diagnoses
    png(sprintf("%s/ptSimilarity.png", path), 3000, 1000)
    heatmap.2(x=simMat,
              dendrogram="both",
              Rowv=TRUE,
              Colv=TRUE,
              ColSideColors = c(brewer.pal(12, "Set3"), brewer.pal(9, "BrBG"))[as.numeric(as.factor(as.character(d)))],
              RowSideColors = c(brewer.pal(12, "Set3"), brewer.pal(9, "BrBG"))[as.numeric(as.factor(as.character(d)))],
              cexRow=0.75,cexCol=1, margins=c(12,12),
              trace="none", key=TRUE,
              col=bluered,
              notecol="black")
    legend("left", legend=unique(sort(as.character(d))), fill=c(brewer.pal(12, "Set3"), brewer.pal(9, "BrBG")), cex=2)
    dev.off()
  }
  return(0)
}


