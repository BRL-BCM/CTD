#' Minimum encoding length (MLE)
#'
#' This function calculates the mininmum encoding length associated with a subset of variables given a background knowledge graph.
#' @param bs - A list of bitstrings associated with a given patient's perturbed variables.
#' @param pvals - The matrix that gives the perturbation strength significance for all variables (columns) for each patient (rows)
#' @param ptID - The row name in data.pvals corresponding to the patient you specifically want encoding information for.
#' @param G - A list of probabilities with list names being the node names of the background graph.
#' @return df - a data.frame object, for every bitstring provided in bs input parameter, a row is returned with the following data:
#'              the patientID; the bitstring evaluated where T denotes a hit and 0 denotes a miss; the subsetSize, or the number of
#'              hits in the bitstring; the individual p-values associated with the variable's perturbations, delimited by '/';
#'              the combined p-value of all variables in the set using Fisher's method; Shannon's entropy, IS.null;
#'              the minimum encoding length IS.alt; and IS.null-IS.alt, the d.score.
#' @export mle.getEncodingLength
#' @keywords minimum length encoding
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # Identify the most significant subset per patient, given the background graph
#' data_mx.pvals = t(apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE)))
#' for (pt in 1:ncol(data_mx)) {
#'     ptID = colnames(data_mx)[pt]
#'     res = mle.getEncodingLength(ptBSbyK[[ptID]], data_mx.pvals, ptID, G)
#'     res = res[order(res[,"d.score"], decreasing=TRUE),]
#'     print(res)
#' }
mle.getEncodingLength = function(bs, pvals, ptID, G) {
  if (is.null(pvals)) {
    results = data.frame(optimalBS=character(), subsetSize=integer(), opt.T=integer(),
                         IS.null=numeric(), IS.alt=numeric(), d.score=numeric(), stringsAsFactors = FALSE)
  } else {
    results = data.frame(patientID=character(), optimalBS=character(), subsetSize=integer(), opt.T=integer(), varPvalue=numeric(),
                         fishers.Info=numeric(), IS.null=numeric(), IS.alt=numeric(), d.score=numeric(), stringsAsFactors = FALSE)
  }
  row = 1
  for (k in 1:length(bs)) {
    optBS = bs[[k]]
    mets.k = names(optBS)[which(optBS==1)]
    found = sum(optBS)-1
    not_found = k-found
    e = (not_found+1)*log2(length(G)) + length(optBS-1)

    optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
    if (!is.null(pvals) && !is.null(ptID)) {
      results[row, "patientID"] = ptID
      results[row, "varPvalue"] = paste(format(pvals[ptID, mets.k], digits=2, width=3), collapse="/")
      results[row, "fishers.Info"] = -log2(stats.fishersMethod(pvals[ptID, mets.k]))
    }
    results[row, "optimalBS"] = optBS.tmp
    results[row, "subsetSize"] = k+1
    results[row, "opt.T"] = found+1
    results[row, "IS.null"] = log2(choose(length(G), k+1))
    results[row, "IS.alt"] = e
    results[row, "d.score"] = log2(choose(length(G), k+1)) - e
    row = row + 1
  }
  return (results)
}

mle.getEncodingLength_nullk = function(bs, G, null.k) {
  results = data.frame(optimalBS=character(), subsetSize=integer(), opt.T=integer(),
                       IS.null=numeric(), IS.alt=numeric(), d.score=numeric(), stringsAsFactors = FALSE)
  optBS = bs[[1]]
  mets.k = names(optBS)[which(optBS==1)]
  found = sum(optBS)
  not_found = null.k-found
  e = (not_found+1)*log2(length(G)) + length(optBS-1)

  optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
  results[1, "optimalBS"] = optBS.tmp
  results[1, "subsetSize"] = null.k
  results[1, "opt.T"] = found
  results[1, "IS.null"] = log2(choose(length(G), null.k))
  results[1, "IS.alt"] = e
  results[1, "d.score"] = log2(choose(length(G), null.k)) - e

  return (results)
}



