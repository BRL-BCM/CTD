#' Minimum encoding length (MLE)
#'
#' This function calculates the mininmum encoding length associated with a subset of variables given a background knowledge graph.
#' @param bs - A list of bitstrings associated with a given patient's perturbed variables.
#' @param pvals - The matrix that gives the perturbation strength significance for all variables (columns) for each patient (rows)
#' @param ptID - The row name in data.pvals corresponding to the patient you specifically want encoding information for.
#' @export mle.getEncodingLength
#' @keywords minimum length encoding
#' @examples
#' mle.getEncodingLength(optBS, data.pvals, ptID)
mle.getEncodingLength = function(bs, pvals, ptID) {
  results = data.frame(patientID=character(), optimalBS=character(), subsetSize=integer(), opt.T=integer(), varPvalue=numeric(),
                       fishers.Info=numeric(), IS.null=numeric(), IS.alt=numeric(), d.score=numeric(), stringsAsFactors = FALSE)
  row = 1
  for (k in 1:length(bs)) {
    optBS = bs[[k]]
    mets.k = names(optBS)[which(optBS==1)]
    p.k = sum(optBS)

    if(p.k==1) {
      e = log2(length(igraphObjectG))
    } else {
      e = log2(length(igraphObjectG)) + log2(p.k) + 1 + (length(optBS)-1)*stats.entropyFunction(optBS[2:length(optBS)])
    }

    optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
    results[row, "patientID"] = ptID
    results[row, "optimalBS"] = optBS.tmp
    results[row, "subsetSize"] = kk
    results[row, "opt.T"] = p.k
    results[row, "varPvalue"] = paste(format(pvals[ptID, mets.k], digits=2, width=3), collapse="/")
    results[row, "fishers.Info"] = -log2(stats.fishersMethod(pvals[ptID,mets.k]))
    results[row, "IS.null"] = log2(choose(length(V(ig)$name), kk))
    results[row, "IS.alt"] = e
    results[row, "d.score"] = log2(choose(length(V(ig)$name), kk)) - e
    row = row + 1
  }

  return (results)
}
