#' Surrogate profiles
#'
#' Fill in a data matrix rank with surrogate profiles., when your data is
#' low n, high p. 
#' @param data - Data matrix with observations (e.g., patient samples)
#'               as columns, features (e.g., metabolites or genes) as rows
#' @param std - The level of variability (standard deviation) around each
#'              observed feature's z-score you want to add to generate the 
#'              surrogate profiles.
#' @param ref_data - Data matrix for healthy control "reference" samples,
#'                   observations (e.g., patient samples) as columns, 
#'                   features (e.g., metabolites or genes) as rows
#' @return data_mx_surr - Data matrix with added surrogate profiles.
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
#'                   which(colnames(data_mx) %in% diag_pts)]
#' data_mx_surr=data.surrogateProfiles(data=diag_data, std=1, ref_data=refs2)
data.surrogateProfiles = function(data, std=1, ref_data=NULL) {
    if (!is.null(ref_data)) {
      ref_data=ref_data[which(rownames(ref_data) %in% rownames(data)),]
      data=data[which(rownames(data) %in% rownames(ref_data)),]
      rpt=ceiling(nrow(data)/ncol(data)/2)
    } else { rpt=ceiling(nrow(data)/ncol(data)) }
    numSurr = ceiling(nrow(data)/2)
    if (numSurr > ncol(data)) {
      # Generate disease surrogates for each unique disease profile.
      data_surr=matrix(NA, nrow=nrow(data), ncol=ncol(data)+ncol(data)*rpt)
      c_col=ncol(data)+1
      for (pt in seq_len(ncol(data))){
        data_surr[,pt]=data[,pt]
        for (rrpt in seq_len(rpt)){
          rr=rnorm(nrow(data),mean=0,sd=std)
          data_surr[,c_col]=as.numeric(data[,pt])+rr
          c_col=c_col+1
        }
      }
      colnames(data_surr)=c(colnames(data),
                            sprintf("disease_surr%d",
                                    seq_len(ncol(data_surr)-ncol(data))))
      rownames(data_surr)=rownames(data)
    } else {
      ind = sample(seq_len(ncol(data)), numSurr)
      data_surr = data[,ind]
    }
    if (!is.null(ref_data)) {
      if (numSurr > ncol(ref_data)) {
        cntl_surr=matrix(NA, nrow=nrow(ref_data), 
                              ncol=ncol(ref_data)+ncol(ref_data)*rpt)
        c_col = ncol(ref_data)+1
        for (pt in seq_len(ncol(ref_data))) {
          cntl_surr[,pt]=ref_data[,pt]
          for (rrpt in seq_len(rpt)) {
            rr=rnorm(nrow(ref_data), mean=0, sd=std)
            cntl_surr[,c_col]=as.numeric(ref_data[,pt])+rr
            c_col=c_col+1
          }
        }
        colnames(cntl_surr)=c(colnames(ref_data), 
                              sprintf("cntl_surr%d",
                                      seq_len(ncol(cntl_surr)-ncol(ref_data))))
        rownames(cntl_surr)=rownames(ref_data)
      } else {
        ind=sample(seq_len(ncol(ref_data)), numSurr)
        cntl_surr=ref_data[,ind]
      }
      data_mx_surr=cbind(data_surr, cntl_surr)
    } else {data_mx_surr = data_surr}
    # Impute metabolites that are NA for all samples in data_mx_surr
    if (!is.null(ref_data)){data_mx_surr=data.imputeData(data_mx_surr,ref_data)
    } else { data_mx_surr = data.imputeData(data_mx_surr, data) }
    # Remove metabolites that do no vary (have a standard deviation of 0)
    var.met = apply(data_mx_surr, 1, sd)
    if (length(which(var.met == 0)) > 0){
      data_mx_surr=data_mx_surr[-which(var.met == 0),]
    }
    rownames(data_mx_surr) = tolower(rownames(data_mx_surr))
    data_mx_surr = apply(data_mx_surr, c(1, 2), as.numeric)
    return(data_mx_surr)
}