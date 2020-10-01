#' Minimum encoding length 
#'
#' This function calculates the mininmum encoding length associated with a
#' subset of variables given a background knowledge graph.
#' @param bs - A list of bitstrings associated with a given patient's
#'             perturbed variables.
#' @param pvals - The matrix that gives the perturbation strength significance
#'                for all variables (columns) for each patient (rows)
#' @param ptID - The row name in data.pvals corresponding to the patient you
#'               specifically want encoding information for.
#' @param G - A list of probabilities with list names being the node names of
#'            the background graph.
#' @return df - a data.frame object, for every bitstring provided in bs input
#' parameter, a row is returned with the following data: the patientID; the 
#' bitstring evaluated where T denotes a hit and 0 denotes a miss; the
#' subsetSize, or the number of hits in the bitstring; the individual p-values
#' associated with the variable's perturbations, delimited by '/'; the combined
#' p-value of all variables in the set using Fisher's method; Shannon's
#' entropy, IS.null; the minimum encoding length IS.alt; and IS.null-IS.alt,
#' the d.score.
#' @export mle.getEncodingLength
#' @keywords minimum length encoding
#' @importFrom gmp chooseZ
#' @examples
#' # Identify the most significantly connected subset for a given patients'
#' # perturbations, given the network G
#' data("Miller2015")
#' data_mx = Miller2015[-c(1,grep("x - ",rownames(Miller2015))),
#'                         grep("IEM", colnames(Miller2015))]
#' data_mx = apply(data_mx, c(1,2), as.numeric)
#' data_pval=t(apply(data_mx,c(1,2),
#'                     function(i)2*pnorm(abs(i),lower.tail=FALSE)))
#' # Choose patient #1's (i.e., IEM_1000's) top 5 perturbed metabolites
#' ptID = colnames(data_mx)[1]
#' S=rownames(data_mx)[order(abs(data_mx[,which(colnames(data_mx)==ptID)]),
#'                             decreasing=TRUE)[seq_len(5)]]
#' # Build a dummy metabolite network for all metabolites in data_mx
#' adj_mat=matrix(0, nrow=nrow(data_mx), ncol=nrow(data_mx))
#' rows=sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' cols=sample(seq_len(ncol(adj_mat)), 0.1*ncol(adj_mat))
#' for (i in rows){for (j in cols){adj_mat[i,j]=rnorm(1,mean=0,sd=1)}}
#' colnames(adj_mat) = rownames(data_mx)
#' rownames(adj_mat) = rownames(data_mx)
#' G = vector("numeric", length=ncol(adj_mat))
#' names(G)=colnames(adj_mat)
#' ranks = list()
#' for (n in seq_len(length(S))) { 
#'     print(sprintf("%d / %d", n, length(S)))
#'     ind = which(names(G)==S[n])
#'     ranks[[n]]=singleNode.getNodeRanksN(ind,G,p1=0.9,thresholdDiff=0.01,
#'                                         adj_mat,S,log2(length(G)),FALSE) 
#' }
#' names(ranks) = S
#' ptBSbyK = mle.getPtBSbyK(S, ranks)
#' res = mle.getEncodingLength(ptBSbyK, data_pval, ptID, G)
#' # Rows with d.scores > 4.32 are of interest. Anything less indicates
#' # no to weak signal.
#' res = res[order(res[,"d.score"], decreasing=TRUE),]
#' print(res)
mle.getEncodingLength = function(bs, pvals, ptID, G) {
    if (is.null(pvals)) {
        results = data.frame(optimalBS=character(),subsetSize=integer(),
                                opt.T=integer(),IS.null=numeric(),
                                IS.alt=numeric(),d.score=numeric(),
                                stringsAsFactors=FALSE)
    } else {
        results = data.frame(patientID=character(),optimalBS=character(),
                                subsetSize=integer(),opt.T=integer(),
                                varPvalue=numeric(),fishers.Info=numeric(),
                                IS.null=numeric(),IS.alt=numeric(),
                                d.score=numeric(),stringsAsFactors=FALSE)
    }
    row=1
    for (k in seq_len(length(bs))) { #Assume k=1 corresponds to subset of size 1
        optBS=bs[[k]]
        mets.k=names(optBS)[which(optBS==1)]
        found=sum(optBS)
        not_found=k-found
        e=(not_found+1)*log2(length(G)) + length(optBS)-1
        optBS.tmp = gsub("1", "T", paste(as.character(optBS), collapse=""))
        if (!is.null(pvals) && !is.null(ptID)) {
            results[row,"patientID"] = ptID
            results[row,"varPvalue"]=paste(format(pvals[ptID, mets.k],
                                                    digits=2,width=3),
                                            collapse="/")
            fishers.pval=-log2(stat.fishersMethod(pvals[ptID,mets.k]))
            results[row,"fishers.Info"]=fishers.pval
        }
        results[row,"optimalBS"]=optBS.tmp
        results[row,"subsetSize"]=k
        results[row,"opt.T"]=found
        results[row,"IS.null"]=log2(chooseZ(length(G), k))
        results[row,"IS.alt"]=e
        results[row,"d.score"]=round(log2(chooseZ(length(G), k)) - e, 3)
        row=row+1
    }
    return (results)
}



