#' Impute missing values as lowest observed value.
#'
#' @param data - Normalized, imputed data. Data matrix with observations as rows, features as columns.
#' @param ref - Reference samples normalized, imputed data.
#' @return imputed.data - Z-transformed data.
#' @importFrom stats runif
#' @export data.imputeData
data.imputeData = function(data, ref) {
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = data[sort(rownames(data)),]
  ref = ref[sort(rownames(ref)),]
  
  imputed.data = data
  for (met in 1:nrow(ref)) {
    rowData = ref[met,]
    if (any(is.na(rowData))) {
      rowData = as.numeric(rowData[-which(is.na(rowData))])
    } else {
      rowData = as.numeric(rowData)
    }
    # Impute using uniform random variable, where a = 0.99*observed minimum, and b = observed minimum
    min_row = min(rowData)
    if (min_row<0) {
      min_row = -1*min_row
      imputed.data[met, is.na(data[met,])] = tryCatch(-1*runif(sum(is.na(data[met,])), min = 0.99*min_row, max= min_row), 
                                                      error = function(e) e, warning=function(w) print(sprintf("%s: met%d", w, met)))
    } else {
      imputed.data[met, is.na(data[met,])] = tryCatch(runif(sum(is.na(data[met,])), min = 0.99*min(rowData), max= min(rowData)), 
                                                      error = function(e) e, warning=function(w) print(sprintf("%s: met%d", w, met)))
    }
  }
  return(imputed.data)
}
