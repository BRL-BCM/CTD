#' Patient similarity using mutual information MLE metric of patients' most modular, perturbed subsets.
#'
#' This function calculates the universal distance between patients, using a mutual information metric, where self-information comes from the minimum encoding length of each patient's encoded modular perturbations in the background knowledge graph.
#' @param p1.optBS - The optimal bitstring associated with patient 1.
#' @param ptID - The identifier associated with patient 1's sample.
#' @param p2.optBS - The optimal bitstring associated with patient 2.
#' @param ptID - The identifier associated with patient 2's sample.
#' @param data_mx - The matrix that gives the perturbation strength (z-scores) for all variables (columns) for each patient (rows).
#' @return patientSim - a similarity matrix, where row and columns are patient identifiers.
#' @export mle.getPtSim
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
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
#'       tmp = mle.getPtSim(ptBSbyK[[ptID]][k], ptID, ptBSbyK[[ptID2]][k], ptID2, data_mx, perms)
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
mle.getPtSim = function(p1.optBS, ptID, p2.optBS, ptID2, data_mx, perms) {
  if (length(p1.optBS) != length(p2.optBS)) {
    print("Make sure subset from patient1 is the same size as subset from patient2.")
    return(0)
  }

  G = vector(mode="list", length=nrow(data_mx))
  names(G) = rownames(data_mx)
  p1.e = mle.getEncodingLength(p1.optBS, NULL, ptID, G)[,"IS.alt"]
  p2.e = mle.getEncodingLength(p2.optBS, NULL, ptID2, G)[,"IS.alt"]

  # Get optimal bitstring for encoding of patient1's union patient2's subsets
  dirSim = vector("numeric", length=length(p1.optBS))
  for (k in 1:length(p1.optBS)) {
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
    dirSim[k] = 1-(length(intersect(p1.dirs_sign, p2.dirs_sign))/(length(p12.dirs)))
  }

  if (is.null(perms)) {
    # Using the "with memory" network walker, get optimal bitstring for encoding of patient1's union patient2's subsets
    p1.sig.nodes = names(sort(abs(data_mx[,ptID]), decreasing = TRUE)[1:kmx])
    p2.sig.nodes = names(sort(abs(data_mx[,ptID2]), decreasing = TRUE)[1:kmx])
    p12.sig.nodes = unique(c(p1.sig.nodes, p2.sig.nodes))
    perms = mle.getPermN_memory(p12.sig.nodes, G)
    p12.e = c()
    p12.optBS = list()
    for (k in 1:length(p1.bs)) {
      p1.bss = p1.bs[[k]]
      p2.bss = p2.bs[[k]]
      p1.sig.nodes = names(sort(abs(data_mx[,ptID]), decreasing = TRUE)[1:k])
      p2.sig.nodes = names(sort(abs(data_mx[,ptID2]), decreasing = TRUE)[1:k])
      p12.sig.nodes = unique(c(p1.sig.nodes, p2.sig.nodes))
      p12.optBS[[k]] = mle.getPtBSbyK(p12.sig.nodes, perms)[[length(p12.sig.nodes)]]
      p12.e[k] = mle.getEncodingLength_nullk(p12.optBS[k], G, length(p12.sig.nodes))[,"IS.alt"]
    }
  } else {
    # Using the "memoryless" network walker, get optimal bitstring for encoding of patient1's union patient2's subsets.
    p12.e = c()
    for (k in 1:length(p1.optBS)) {
      p1.sig.nodes = names(sort(abs(data_mx[,ptID]), decreasing = TRUE))[1:k]
      p2.sig.nodes = names(sort(abs(data_mx[,ptID2]), decreasing = TRUE))[1:k]
      p12.sig.nodes = unique(c(p1.sig.nodes, p2.sig.nodes))
      p12.optBS = mle.getPtBSbyK_memoryless(p12.sig.nodes, perms)[[length(p12.sig.nodes)]]
      p12.e[k] = mle.getEncodingLength_nullk(list(p12.optBS), G, length(p12.sig.nodes))[,"IS.alt"]
    }
  }

  # Normalized Compression Distance, Percent Mutual Information
  return (list(p1.e=p1.e, p2.e=p2.e, p12.e=p12.e,
               NCD=(p12.e-apply(cbind(p1.e, p2.e), 1, min))/apply(cbind(p1.e, p2.e), 1, max),
               mutualInfoPer=1-((p1.e+p2.e-p12.e)/p12.e),
               dirSim=dirSim))
}

