#' Minimum encoding length (MLE)
#'
#' This function calculates the mininmum encoding length associated with a subset of variables given a background knowledge graph.
#' @param bs - A list of bitstrings associated with a given patient's perturbed variables.
#' @param G - A list of probabilities with list names being the node names of the background graph.
#' @return df - a data.frame object, for every bitstring provided in bs input parameter, a row is returned with the following data:
#'              the patientID; the bitstring evaluated where T denotes a hit and 0 denotes a miss; the subsetSize, or the number of
#'              hits in the bitstring; the individual p-values associated with the variable's perturbations, delimited by '/';
#'              the combined p-value of all variables in the set using Fisher's method; Shannon's entropy, IS.null;
#'              the minimum encoding length IS.alt; and IS.null-IS.alt, the d.score.
#' @export kraft.getEncodingLength
#' @keywords minimum length encoding
kraft.getEncodingLength_memoryless = function(bs, G) {
  results = data.frame(optimalBS=character(), subsetSize=integer(), opt.T=integer(),
                       IS.null=numeric(), IS.alt=numeric(), d.score=numeric(), stringsAsFactors = FALSE)
  row = 1
  for (k in 1:length(bs)) {
    optBS = bs[[k]]
    mets.k = names(optBS)[which(optBS==1)]
    found = sum(optBS)
    not_found = k+1-found

    e = log2(length(G)) + log2(choose(length(G), not_found)) + stats.iteratedLog2(found-1)
    if (length(optBS)>1) {
      e = e + stats.iteratedLog2(length(optBS)-1) + (length(optBS)-2)*stats.entropyFunction(optBS[2:(length(optBS)-1)])
    }

    optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
    results[row, "patientID"] = ptID
    results[row, "optimalBS"] = optBS.tmp
    results[row, "subsetSize"] = k+1
    results[row, "opt.T"] = found
    results[row, "varPvalue"] = paste(format(pvals[ptID, mets.k], digits=2, width=3), collapse="/")
    results[row, "fishers.Info"] = -log2(stats.fishersMethod(pvals[ptID, mets.k]))
    results[row, "IS.null"] = log2(choose(length(G), k+1))
    results[row, "IS.alt"] = e
    results[row, "d.score"] = log2(choose(length(G), k+1)) - e
    row = row + 1
  }

  return (results)
}


mle.getEncodingLength_memory = function(bs, pvals, ptID, G) {
  results = data.frame(patientID=character(), optimalBS=character(), subsetSize=integer(), opt.T=integer(), varPvalue=numeric(),
                       fishers.Info=numeric(), IS.null=numeric(), IS.alt=numeric(), d.score=numeric(), stringsAsFactors = FALSE)
  row = 1
  for (k in 1:length(bs)) {
    optBS = bs[[k]]
    mets.k = names(optBS)[which(optBS==1)]
    found = sum(optBS)-1
    not_found = k-found

    e = log2(length(G)) + log2(choose(length(G), not_found))
    hits = which(optBS==1)[-1]
    if (length(hits)>0) {
      prevHit = 1
      onesIsland=1
      for (h in 1:length(hits)) {
        numZeros = hits[h]-prevHit
        if (numZeros>1) {
          e = e + stats.iteratedLog2(numZeros)
          if (onesIsland>1) {
            e = e + stats.iteratedLog2(onesIsland)
          } else {
            e = e + 1
          }
          onesIsland = 1
        } else {
          e = e + 1
          onesIsland = onesIsland + 1
        }
        prevHit = hits[h]
      }
    }

    optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
    results[row, "patientID"] = ptID
    results[row, "optimalBS"] = optBS.tmp
    results[row, "subsetSize"] = k+1
    results[row, "opt.T"] = found
    results[row, "varPvalue"] = paste(format(pvals[ptID, mets.k], digits=2, width=3), collapse="/")
    results[row, "fishers.Info"] = -log2(stats.fishersMethod(pvals[ptID, mets.k]))
    results[row, "IS.null"] = log2(choose(length(G), k+1))
    results[row, "IS.alt"] = e
    results[row, "d.score"] = log2(choose(length(G), k+1)) - e
    row = row + 1
  }

  return (results)
}
