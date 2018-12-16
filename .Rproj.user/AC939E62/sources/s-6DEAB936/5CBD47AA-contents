#' Generate heatmap plot of patient similarity matrix.
#'
#' This function plots a heatmap of a patient similarity matrix.
#' @param simMat - The patient similarity matrix.
#' @param path - The filepath to a directory in which you want to store the .png file.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of simMat.
#' @export plot.hmSim
#' @examples
#' # Read in any network via its adjacency matrix
#' tmp = as.matrix(read.table("adjacency_matrix.txt", sep="\t", header=TRUE))
#' colnames(tmp) = rownames(tmp)
#' ig = graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
#' V(ig)$name = tolower(V(ig)$name)
#' adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))  # Must have this declared as a GLOBAL variable!!!!!
#' # Set other tuning parameters
#' p0=0.1  # 10% of probability distributed uniformly
#' p1=0.9  # 90% of probability diffused based on edge weights in networks
#' G = vector(mode="list", length=length(V(ig)$name))
#' names(G) = names(V(ig)$name)
#' # Get node permutations for graph
#' perms = list()
#' for (n in 1:length(G)) {
#'     print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
#'     perms[[names(G)[n]]] = mle.getPermN(n, G)
#' }
#' # Decide what the largest subset size you will consider will be
#' kmx = 20
#' # Load your patient data (p features as rows x n observations as columns)
#' # data_mx = read.table("/your/own/data.txt", sep="\t", header=TRUE)
#' data(testData)
#' data_mx = testData
#' # Get bitstrings associated with each patient's top kmx variable subsets
#' ptBSbyK = list()
#' for (pt in 1:ncol(data_mx)) {
#'     ptID = colnames(data_mx)[pt]
#'     ptBSbyK[[ptID]] = mle.getPtBSbyK(data_mx, ptID, perms, kmx)
#' }
#' # Get patient distances
#' patientSim = matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx))
#' rownames(patientSim) = colnames(data_mx)
#' colnames(patientSim) = colnames(data_mx)
#' for (pt in 1:ncol(data_mx)) {
#'     ptID = colnames(data_mx)[pt]
#'     for (pt2 in 1:ncol(data_mx)) {
#'         ptID2 = colnames(data_mx)[pt2]
#'         patientSim[ptID, ptID2] = mle.getPatientSimilarity(ptBSbyK[ptID], ptID, ptBSbyK[ptID2], data_mx)
#'     }
#' }
#' # if you have diagnostic labels associated with the colnames(data_mx), send them using diagnoses parameter
#' diagnoses = colnames(data_mx)
#' diagnoses[1:50] = "diseased"
#' diagnoses[51:100] = "neg_control"
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


