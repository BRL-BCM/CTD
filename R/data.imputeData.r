#' Impute missing values
#' 
#' Impute missing values as lowest observed value in a reference population
#' @param data - Normalized data with some missingness. Data matrix with 
#'               features as rows, samples as columns.
#' @param ref - Reference sample data with features as rows, samples as
#'              columns. Can include some missingness.
#' @return imputed.data - Imputed data.
#' @importFrom stats runif
#' @export data.imputeData
#' @examples
#' data(Thistlethwaite2020)
#' data_mx = Thistlethwaite2020
#' # Data with missing values
#' dt_w_missing_vals = data_mx[1:25,-seq_len(8)]
#' # Reference data can also have missing values
#' ref_data = data_mx[1:25,grep("EDTA-REF", colnames(data_mx))]
#' fil.rate = apply(ref_data, 1, function(i) sum(is.na(i))/length(i))
#' # Can only impute data that are found in reference samples
#' dt_w_missing_vals = dt_w_missing_vals[which(fil.rate<1.0),]
#' ref_data = ref_data[which(fil.rate<1.0),]
#' imputed.data = data.imputeData(dt_w_missing_vals, ref_data)
#' print(any(is.na(imputed.data)))
data.imputeData = function(data, ref) {
    # Remove metabolites that are NA for all data/ref samples
    rmThese.ref = c()
    for(r in seq_len(nrow(ref))){
        if(all(is.na(as.numeric(ref[rownames(ref)[r],])))){
            rmThese.ref = c(rmThese.ref, r)}}
    if(length(rmThese.ref)>0){ref=ref[-rmThese.ref,]}
    rmThese.data = c()
    for (r in seq_len(nrow(data))) {
        if (all(is.na(as.numeric(data[rownames(data)[r],])))){
            rmThese.data=c(rmThese.data, r)}}
    if(length(rmThese.data)>0){data=data[-rmThese.data,]}
    # Match metabolites between data and ref matrices
    data = data[which(rownames(data) %in% rownames(ref)),]
    ref = ref[which(rownames(ref) %in% rownames(data)),]
    data = data[sort(rownames(data)),]
    ref = ref[sort(rownames(ref)),]
    imputed.data = data
    for (met in seq_len(nrow(ref))) {
        rowData = ref[met,]
        if (any(is.na(rowData))) {
            rowData = as.numeric(rowData[-which(is.na(rowData))])
        } else {rowData = as.numeric(rowData)}
        # Impute using uniform random variable, where 
        # a = 0.99*observed minimum, and b = observed minimum
        min_row = min(rowData)
        cols = which(is.na(data[met,]))
        if (min_row<0) {
            min_row = -1*min_row
            i_val = -1
        } else {i_val = 1}
        imputed.data[met,cols]=tryCatch(i_val*runif(length(cols),
                                                    min = 0.99*min_row,
                                                    max= min_row),
                                        error = function(e) e,
                                        warning=function(w) w)
    }
    return(imputed.data)
}
