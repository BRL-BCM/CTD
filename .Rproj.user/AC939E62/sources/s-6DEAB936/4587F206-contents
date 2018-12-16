#' Minimum encoding length (MLE)
#'
#' This function calculates the mininmum encoding length associated with a subset of variables given a background knowledge graph.
#' @param bs - A list of bitstrings associated with a given patient's perturbed variables.
#' @param pvals - The matrix that gives the perturbation strength significance for all variables (columns) for each patient (rows)
#' @param ptID - The row name in data.pvals corresponding to the patient you specifically want encoding information for.
#' @export mle.getEncodingLength
#' @keywords minimum length encoding
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
#' # Identify the most significant subset per patient, given the background graph
#' data_mx.pvals = apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
#' for (pt in 1:ncol(data_mx)) {
#'     ptID = colnames(data_mx)[pt]
#'     res = mle.getEncodingLength(ptBSbyK[[ptID]], data_mx.pvals, ptID)
#'     res = res[which.max(res[,"d.score"]),]
#'     print(res)
#' }
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
      e = log2(length(igraphObjectG)) + log2(p.k+1) + (length(optBS)-1)*stats.entropyFunction(optBS[2:length(optBS)])
    }

    optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
    results[row, "patientID"] = ptID
    results[row, "optimalBS"] = optBS.tmp
    results[row, "subsetSize"] = k
    results[row, "opt.T"] = p.k
    results[row, "varPvalue"] = paste(format(pvals[ptID, mets.k], digits=2, width=3), collapse="/")
    results[row, "fishers.Info"] = -log2(stats.fishersMethod(pvals[ptID, mets.k]))
    results[row, "IS.null"] = log2(choose(length(igraphObjectG), k))
    results[row, "IS.alt"] = e
    results[row, "d.score"] = log2(choose(length(igraphObjectG), k)) - e
    row = row + 1
  }

  return (results)
}
