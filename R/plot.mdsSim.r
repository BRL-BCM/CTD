#' View patient clusters using multi-dimensional scaling.
#'
#' This function plots the provided patient similarity matrix in a lower dimensional space
#' using multi-dimensional scaling, which is well suited for similarity metrics.
#' @param simMat - The patient similarity matrix.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of simMat.
#' @param k - The number of dimension you want to plot your data using multi-dimensional scaling.
#' @param diag - The diagnosis associated with positive controls in your data.
#' @export plot.mdsSim
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
#' p = plot.mdsSim(patientSim, diagnoses, k=2, diag="diseased")
#' p
#' p = plot.mdsSim(patientSim, diagnoses, k=3, diag="diseased")
#' p
plot.mdsSim = function(simMat, diagnoses, k, diag) {
  if (!(k %in% c(2,3))) {
    print("K must be either 2-dimensions or 3-dimensions.")
    return(0)
  }

  if (is.null(diagnoses)) {
    print("To view patient clusters, please provide clinical labels, corresponding to each patient in the provided similarity matrix.")
    return(0)
  }

  fitSim = cmdscale(patientSim, eig=FALSE, k=k)
  x = round(fitSim[,1], 2)
  y = round(fitSim[,2], 2)
  if (k==3) {
    z = round(fitSim[,3], 2)
    df = data.frame(x=x, y=y, z=z, color=diagnoses, label=colnames(patientSim))
    p = plot_ly(df, x=~x, y=~y, z=~z, color=~color, text=~label, marker = list(size = 20))
  } else {
    df = data.frame(x=x, y=y, color=diagnoses, label=colnames(patientSim))
    p = plot_ly(df, x=~x, y=~y, color=~color, text=~label, marker = list(size = 20))
  }

  # Do K-means here
  #dis.centroid = apply(patientSim[which(diagnoses==diag),], 2, mean)
  #negCntl.centroid = apply(patientSim[which(diagnoses=="negCntl"),], 2, mean)
  #kmns = kmeans(patientSim, centers = rbind(dis.centroid, negCntl.centroid))
  #dcntd = diagnoses[which(colnames(patientSim) %in% names(kmns$cluster[which(kmns$cluster==1)]))]
  #ncntd = diagnoses[which(colnames(patientSim) %in% names(kmns$cluster[which(kmns$cluster==2)]))]

  #rownames(df) = df[,"label"]
  #min(df[dcntd,"x"])
  #max(df[dcntd,"x"])
  #min(df[dcntd,"y"])
  #max(df[dcntd,"y"])

  #min(df[ncntd,"x"])
  #max(df[ncntd,"x"])
  #min(df[ncntd,"y"])
  #max(df[ncntd,"y"])

  # https://plot.ly/r/shapes/
  # Circles are centered around ((x0+x1)/2, (y0+y1)/2))
  #p = layout(p, title = sprintf("Kmeans Clustering of Patient Distances from MolPhenoMatch"),
  #           shapes = list(
  #             list(type = 'circle',
  #                  xref = 'x', x0 = .2, x1 = .7,
  #                  yref = 'y', y0 = 20, y1 = 3000,
  #                  fillcolor = 'rgb(50, 20, 90)', line = list(color = 'rgb(50, 20, 90)'),
  #                  opacity = 0.2),
  #             list(type = 'circle',
  #                  xref = 'x', x0 = .75, x1 = 1.5,
  #                  yref = 'y', y0 = 2500, y1 = 7500,
  #                  fillcolor = 'rgb(30, 100, 120)', line = list(color = 'rgb(30, 100, 120)'),
  #                  opacity = 0.2)))

  return(p)
}
