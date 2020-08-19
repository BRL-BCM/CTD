#' Combine datasets from different source files.
#'
#' @param curr_data - Current data matrix
#' @param more_data - Data matrix you want to combine with curr_data.
#' @return combined.data - Combined data matrix.
#' @export data.combineData
#' @examples 
#' # Row names and column names are required for both curr_data and data input matrices.
#' curr_data = matrix(rnorm(500), ncol=100)
#' rownames(curr_data) = sprintf("Feature%d", sample(seq_len(20), nrow(curr_data), replace = FALSE))
#' colnames(curr_data) = sprintf("Sample%d", seq_len(ncol(curr_data)))
#' more_data = matrix(rnorm(500), ncol=100)
#' rownames(more_data) = sprintf("Feature%d", sample(seq_len(20), nrow(curr_data), replace = FALSE))
#' colnames(more_data) = sprintf("Sample%d", seq_len(ncol(curr_data)))
#' combined.data = data.combineData(curr_data, more_data)
data.combineData = function(curr_data, more_data) {
  curr_data = as.matrix(curr_data)
  more_data = as.matrix(more_data)
  
  unionMets = unique(c(rownames(curr_data), rownames(more_data)))
  combined.data = matrix(0, nrow=length(unionMets), ncol=ncol(curr_data)+ncol(more_data))
  rownames(combined.data) = unionMets
  colnames(combined.data) = c(colnames(curr_data), colnames(more_data))
  for (r in seq_len(length(unionMets))) {
    if (unionMets[r] %in% rownames(curr_data)) {
      combined.data[r,colnames(curr_data)] = as.numeric(curr_data[unionMets[r], ])
    } else {
      combined.data[r,colnames(curr_data)] = rep(NA, ncol(curr_data))
    }
    
    if (unionMets[r] %in% rownames(more_data)) {
      combined.data[r,colnames(more_data)] = as.numeric(more_data[unionMets[r], ])
    } else {
      combined.data[r,colnames(more_data)] = rep(NA, ncol(more_data))
    }
  }
  
  return(combined.data)
}
