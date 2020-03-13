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
#' # Get patient distances
#' data_mx.pvals = apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
#' res = list()
#' tt = list(ncd=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
#' rownames(t$ncd) = colnames(data_mx)
#' colnames(t$ncd) = colnames(data_mx)
#' for (i in 1:kmx) {
#'   res[[i]] = tt
#' }
#' for (pt in 1:ncol(data_mx)) {
#'   print(pt)
#'   ptID = colnames(data_mx)[pt]
#'   for (pt2 in pt:ncol(data_mx)) {
#'     ptID2 = colnames(data_mx)[pt2]
#'     for (k in 1:kmx) {
#'       tmp = mle.getPtSim(ptBSbyK[[ptID]][k], ptID, ptBSbyK[[ptID2]][k], ptID2, data_mx, ranks)
#'       res[[k]]$ncd[ptID, ptID2] = tmp$NCD
#'       res[[k]]$ncd[ptID2, ptID] = tmp$NCD
#'     }
#'   }
#' }
mle.getPtSim = function(p1.optBS, ptID, p2.optBS, ptID2, data_mx, ranks) {
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
    p1.sig.nodes = rownames(data_mx)[order(abs(data_mx[, ptID]), decreasing = TRUE)][1:k]
    p2.sig.nodes = rownames(data_mx)[order(abs(data_mx[, ptID2]), decreasing = TRUE)][1:k]
    p1.dirs = data_mx[p1.sig.nodes, ptID]
    p1.dirs[which(!(p1.dirs>0))] = 0
    p1.dirs[which(p1.dirs>0)] = 1
    p2.dirs = data_mx[p2.sig.nodes, ptID2]
    p2.dirs[which(!(p2.dirs>0))] = 0
    p2.dirs[which(p2.dirs>0)] = 1
    p1.sig.nodes = sprintf("%s%d", p1.sig.nodes, p1.dirs)
    p2.sig.nodes = sprintf("%s%d", p2.sig.nodes, p2.dirs)
    dirSim[k] = 1 - (length(intersect(p1.sig.nodes, p2.sig.nodes))/length(union(p1.sig.nodes, p2.sig.nodes)))
  }

  if (is.null(ranks)) {
    # Using the single-node network walker, get optimal bitstring for encoding of patient1's union patient2's subsets
    p1.sig.nodes = sapply(names(sort(abs(data_mx[,ptID]), decreasing = TRUE)[1:length(p1.optBS)]), trimws)
    p2.sig.nodes = sapply(names(sort(abs(data_mx[,ptID2]), decreasing = TRUE)[1:length(p2.optBS)]), trimws)
    p12.sig.nodes = as.character(sapply(unique(c(p1.sig.nodes, p2.sig.nodes)), trimws))
    ranks = list()
    for (i in 1:length(p12.sig.nodes)) {
      ind = which(names(G)==p12.sig.nodes[i])
      ranks[[i]] = singleNode.getNodeRanksN(n=ind, G=G, S=p12.sig.nodes, num.misses=log2(length(G)))
    }
    p12.e = c()
    for (k in 1:length(p1.optBS)) {
      p1.sig.nodes_cpy = p1.sig.nodes
      p2.sig.nodes_cpy = p2.sig.nodes

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
      p12.sig.nodes_k = sapply(unique(c(p1.sig.nodes_k, p2.sig.nodes_k)), trimws)
      p12.optBS = singleNode.getPtBSbyK(p12.sig.nodes_k, ranks)
      res = mle.getEncodingLength(p12.optBS, NULL, NULL, G)
      p12.e[k] = res[which.max(res[,"d.score"]), "IS.alt"] + log2(choose(length(G), 1))*(length(p12.sig.nodes_k)-which.max(res[,"d.score"]))
    }
  } else {
    # Using predefined node ranks, get optimal bitstring for encoding of patient1's union patient2's subsets.
    p1.sig.nodes = names(sort(abs(data_mx[,ptID]), decreasing = TRUE)[1:length(p1.optBS)])
    p2.sig.nodes = names(sort(abs(data_mx[,ptID2]), decreasing = TRUE)[1:length(p2.optBS)])
    p12.sig.nodes = unique(c(p1.sig.nodes, p2.sig.nodes))
    
    p12.e = c()
    for (k in 1:length(p1.optBS)) {
      p1.sig.nodes_cpy = p1.sig.nodes
      p2.sig.nodes_cpy = p2.sig.nodes
      
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
      p12.sig.nodes_k = sapply(unique(c(p1.sig.nodes_k, p2.sig.nodes_k)), trimws)
      p12.optBS = singleNode.getPtBSbyK(p12.sig.nodes_k, ranks)
      res = mle.getEncodingLength(p12.optBS, NULL, NULL, G)
      p12.e[k] = res[which.max(res[,"d.score"]), "IS.alt"] + log2(choose(length(G), 1))*(length(p12.sig.nodes_k)-which.max(res[,"d.score"]))
    }
  }

  # Normalized Compression Distance, Percent Mutual Information, Jaccard Set Similarity (w/ Directionality)
  return (list(p1.e=p1.e, p2.e=p2.e, p12.e=p12.e,
               NCD=(p12.e-apply(cbind(p1.e, p2.e), 1, min))/apply(cbind(p1.e, p2.e), 1, max),
               mutualInfoPer=1-((p1.e+p2.e-p12.e)/p12.e),
               dirSim=dirSim))
}

