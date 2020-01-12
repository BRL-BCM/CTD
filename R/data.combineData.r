#' Combine datasets from different source files.
#'
#' @param data - Normalized, imputed data. Data matrix with observations as rows, features as columns.
#' @param ref - Reference samples normalized, imputed data.
#' @return combined.data - Z-transformed data.
#' @export data.combineData
data.combineData = function(ref, data) {
  ref = as.matrix(ref)
  data = as.matrix(data)
  
  unionMets = unique(c(rownames(ref), rownames(data)))
  combined.data = matrix(0, nrow=length(unionMets), ncol=ncol(ref)+ncol(data))
  rownames(combined.data) = unionMets
  colnames(combined.data) = c(colnames(ref), colnames(data))
  for (r in 1:length(unionMets)) {
    if (unionMets[r] %in% rownames(ref)) {
      combined.data[r,colnames(ref)] = as.numeric(ref[unionMets[r], ])
    } else {
      combined.data[r,colnames(ref)] = rep(NA, ncol(ref))
    }
    
    if (unionMets[r] %in% rownames(data)) {
      combined.data[r,colnames(data)] = as.numeric(data[unionMets[r], ])
    } else {
      combined.data[r,colnames(data)] = rep(NA, ncol(data))
    }
  }
  
  return(combined.data)
}
