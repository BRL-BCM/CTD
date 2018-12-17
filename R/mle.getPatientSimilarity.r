#' Patient similarity using mutual information MLE metric of patients' most modular, perturbed subsets.
#'
#' This function calculates the universal distance between patients, using a mutual information metric, where self-information comes from the minimum encoding length of each patient's encoded modular perturbations in the background knowledge graph.
#' @param p1.optBS - The optimal bitstring associated with patient 1.
#' @param ptID - The identifier associated with patient 1's sample.
#' @param p2.optBS - The optimal bitstring associated with patient 2.
#' @param ptID - The identifier associated with patient 2's sample.
#' @param data_mx - The matrix that gives the perturbation strength (z-scores) for all variables (columns) for each patient (rows).
#' @export mle.getPatientSimilarity
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
mle.getPatientSimilarity = function(p1.optBS, ptID, p2.optBS, ptID2, data_mx, perms) {
  p1.bs = unlist(p1.optBS[[1]])
  p1.sig.nodes = names(p1.bs[which(p1.bs==1)])
  p2.bs = unlist(p2.optBS[[1]])
  p2.sig.nodes = names(p2.bs[which(p2.bs==1)])

  # Get optimal bitstring for encoding of patient1's union patient2's subsets
  p12.sig.nodes = unique(c(p1.sig.nodes, p2.sig.nodes))
  bitStrings.pt = list()
  for (sn in 1:length(p12.sig.nodes)) {
    startNode = p12.sig.nodes[sn]
    current_node_set = perms[[startNode]]
    bitStrings.pt[[sn]] = as.numeric(current_node_set %in% p12.sig.nodes)
    names(bitStrings.pt[[sn]]) = current_node_set
  }
  tmp = unlist(lapply(bitStrings.pt, function(i) sum(i)))
  bitStrings.pt = bitStrings.pt[which(tmp==max(tmp))]
  p12.optBS = bitStrings.pt[[which.min(unlist(lapply(bitStrings.pt, function(i) sum(which(i==1)))))]]
  ind = which(sapply(p12.optBS, function(i) i==1))
  p12.optBS = p12.optBS[1:ind[length(ind)]]

  p1.e = log2(length(perms)) + log2(length(p1.sig.nodes)+1) + (length(p1.bs)-1)*stats.entropyFunction(p1.bs[2:length(p1.bs)])
  p2.e = log2(length(perms)) +  log2(length(p2.sig.nodes)+1) + (length(p2.bs)-1)*stats.entropyFunction(p2.bs[2:length(p2.bs)])
  p12.e = log2(length(perms)) +  log2(length(p12.sig.nodes)+1) + (length(p12.optBS)-1)*stats.entropyFunction(p12.optBS[2:length(p12.optBS)])

  p1.dirs = data_mx[p1.sig.nodes,ptID]
  p1.dirs[p1.dirs>0] = 1
  p1.dirs[p1.dirs<0] = 0
  names(p1.dirs) = p1.sig.nodes
  p2.dirs = data_mx[p2.sig.nodes,ptID2]
  p2.dirs[p2.dirs>0] = 1
  p2.dirs[p2.dirs<0] = 0
  names(p2.dirs) = p2.sig.nodes

  p1.dirs_sign = sprintf("%s%d", names(p1.dirs), p1.dirs)
  p2.dirs_sign = sprintf("%s%d", names(p2.dirs), p2.dirs)
  p12.dirs = c(p1.dirs, p2.dirs)
  ind = which(duplicated(sprintf("%s%d", names(p12.dirs), p12.dirs)))
  if (length(ind)>0) {
    p12.dirs = p12.dirs[-ind]
  }
  dirSim = 1-(length(intersect(p1.dirs_sign, p2.dirs_sign))/(length(p12.dirs)))
  ncdSim = max(c(p12.e-p1.e, p12.e-p2.e))/max(c(p1.e, p2.e))
  mutSim = 1-((p1.e+p2.e-p12.e)/p12.e)

  # Normalized Compression Distance, Percent Mutual Information
  return (list(p1.e=p1.e, p2.e=p2.e, p12.e=p12.e, p1.k=length(p1.sig.nodes), p2.k=length(p2.sig.nodes), p12.k=length(p12.sig.nodes),
               NCD=ncdSim, mutualInfoPer=mutSim, dirSim=dirSim))
}
