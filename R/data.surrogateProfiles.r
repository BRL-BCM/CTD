#' Generate surrogate profiles
#'
#' Fill in a data matrix rank with surrogate profiles., when your data is
#' low n, high p. 
#' @param data - Data matrix with observations (e.g., patient samples)
#'               as columns, features (e.g., metabolites or genes) as rows
#' @param std - The level of variability (standard deviation) around each
#'              observed feature's z-score you want to add to generate the 
#'              surrogate profiles
#' @param ref_data - Data matrix for healthy control "reference" samples,
#'                   observations (e.g., patient samples) as columns, 
#'                   features (e.g., metabolites or genes) as rows
#' @return data_mx_surr - Data matrix with added surrogate profiles
#' @importFrom stats na.omit rnorm sd
#' @export data.surrogateProfiles
#' @examples
#' data("Miller2015")
#' data_mx=Miller2015[-1,grep("IEM_", colnames(Miller2015))]
#' data_mx=apply(data_mx, c(1,2), as.numeric)
#' diags=unlist(Miller2015["diagnosis",grep("IEM_", colnames(Miller2015))])
#' refs=data_mx[,which(diags=="No biochemical genetic diagnosis")]
#' ref_fill=as.numeric(Miller2015$`Times identifed in all 200 samples`[-1])/200
#' refs2=refs[which(ref_fill>0.8),]
#' diag_pts=names(diags[which(diags==unique(diags)[1])])
#' diag_data=data_mx[which(rownames(data_mx) %in% rownames(refs2)), 
#'                     which(colnames(data_mx) %in% diag_pts)]
#' data_mx_surr=data.surrogateProfiles(data=diag_data, std=1, ref_data=refs2)
data.surrogateProfiles = function(data, std=1, ref_data) {
    ref_data=ref_data[which(rownames(ref_data) %in% rownames(data)),]
    data=data[which(rownames(data) %in% rownames(ref_data)),]
    rpt=ceiling(nrow(data)/ncol(data)/2)
    numSurr = ceiling(nrow(data)/2)
    if (numSurr>ncol(data)) { # Generate disease surrogates
        d_surr=matrix(NA, nrow=nrow(data),ncol=ncol(data)+ncol(data)*rpt)
        c_col=ncol(data)+1
        for (pt in seq_len(ncol(data))){
            d_surr[,pt]=data[,pt]
            for (rrpt in seq_len(rpt)){
                rr=rnorm(nrow(data),mean=0,sd=std)
                d_surr[,c_col]=as.numeric(data[,pt])+rr
                c_col=c_col+1}}
        colnames(d_surr)=c(colnames(data),
                            sprintf("disease_surr%d",
                                    seq_len(ncol(d_surr)-ncol(data))))
        rownames(d_surr)=rownames(data)} else {d_surr=data}
    if (numSurr>ncol(ref_data)) { # Generate control surrogates
        c_surr=matrix(NA, nrow=nrow(ref_data),
                      ncol=ncol(ref_data)+ncol(ref_data)*rpt)
        c_col = ncol(ref_data)+1
        for (pt in seq_len(ncol(ref_data))) {
            c_surr[,pt]=ref_data[,pt]
            for (rrpt in seq_len(rpt)) {
                rr=rnorm(nrow(ref_data), mean=0, sd=std)
                c_surr[,c_col]=as.numeric(ref_data[,pt])+rr
                c_col=c_col+1}}
        colnames(c_surr)=c(colnames(ref_data),
                            sprintf("control_surr%d",
                                    seq_len(ncol(c_surr)-ncol(ref_data))))
        rownames(c_surr)=rownames(ref_data)} else {c_surr=ref_data}
    d_surr=d_surr[,c(seq_len(ncol(data)),
                    sample(seq(ncol(data)+1,ncol(d_surr)),
                            numSurr-ncol(data)))]
    c_surr=c_surr[,c(seq_len(ncol(ref_data)),
                    sample(seq(ncol(ref_data)+1,ncol(c_surr)),
                            numSurr-ncol(data)))]
    data_mx_surr=cbind(d_surr, c_surr)
    # Impute metabolites that are NA
    if (!is.null(ref_data)){data_mx_surr=data.imputeData(data_mx_surr,ref_data)
    } else { data_mx_surr = data.imputeData(data_mx_surr, data) }
    var.met = apply(data_mx_surr, 1, sd) # Remove metabolites that do no vary
    if (length(which(var.met == 0)) > 0){
        data_mx_surr=data_mx_surr[-which(var.met == 0),]}
    rownames(data_mx_surr) = tolower(rownames(data_mx_surr))
    data_mx_surr = apply(data_mx_surr, c(1, 2), as.numeric)
    return(data_mx_surr)
}