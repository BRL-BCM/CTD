#' Z-transform available data
#' 
#' The z-transform is meant to work with normalized,
#' imputed metabolomics data
#' @param data - Normalized, imputed data. Data matrix includes 
#'               features as rows, samples as columns.
#' @param ref - Normalized, imputed reference sample data. Data
#'              includes features as rows, samples as columns.
#' @return zscored.data - Z-transformed data.
#' @importFrom stats quantile qqnorm lm
#' @export data.zscoreData
#' @examples 
#' dis_data = matrix(rexp(500), ncol=100)
#' rownames(dis_data)=sprintf("Feature%d",seq_len(nrow(dis_data)))
#' colnames(dis_data)=sprintf("Sample%d",seq_len(ncol(dis_data)))
#' ref_data = matrix(rexp(500), ncol=100)
#' rownames(ref_data)=sprintf("Feature%d",seq_len(nrow(ref_data)))
#' colnames(ref_data)=sprintf("Sample%d",seq_len(ncol(ref_data)))
#' zscored.data=data.zscoreData(dis_data,ref_data)
data.zscoreData = function(data, ref) {
    print("zscoreData() called.")
    # Only metabolites found in the reference population can be z-scored
    data = data[which(rownames(data) %in% rownames(ref)),]
    ref = ref[which(rownames(ref) %in% rownames(data)),]
    data = data[sort(rownames(data)),]
    ref = ref[sort(rownames(ref)),]
    # Log transform data
    data = log(data)
    ref = log(data.matrix(ref))
    zscored.data = data
    for (met in seq_len(nrow(data))) {
        met_data = as.numeric(ref[met,])
        rmSamples=unique(c(which(is.na(met_data)),which(is.infinite(met_data))))
        if (length(rmSamples)>0) {x=met_data[-rmSamples]} else {x=met_data}
        if (!all(is.na(x))) {
            if (length(x[intersect(which(x>quantile(x, 0.025)), 
                                    which(x<quantile(x, .975)))])>3) {
                x = x[intersect(which(x>quantile(x, 0.025)),
                                which(x<quantile(x, .975)))]
            }
            d = qqnorm(x, plot.it = FALSE);
            x = as.numeric(d$y)
            z = as.numeric(d$x)
            df = data.frame(x=x,z=z)
            t = lm(x~z, data=df)
            mn.est = as.numeric(t$coefficients[1])
            sd.est = as.numeric(t$coefficients[2])
            rm(d,x,z,df,t)
            zscored.data[met,] = (data[met, ]-mn.est)/sd.est
        }
    }
    return(zscored.data)
}
