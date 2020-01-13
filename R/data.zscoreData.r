#' Z-transform normalized, imputed metabolomics data
#'
#' @param data - Normalized, imputed data. Data matrix with observations as rows, features as columns.
#' @param ref - Reference samples normalized, imputed data.
#' @return zscore.data - Z-transformed data.
#' @export data.zscoreData
data.zscoreData = function(data, ref) {
  print("zscoreData() called.")
  
  # Only metabolites that also occur in the reference population can be z-scored
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = data[sort(rownames(data)),]
  ref = ref[sort(rownames(ref)),]
  
  # Log transform data
  data = log(data)
  ref = log(data.matrix(ref))
  
  zscore.data = data
  for (met in 1:nrow(data)) {
    met_data = as.numeric(ref[met,])
    rmSamples = unique(c(which(is.na(met_data)), which(is.infinite(met_data))))
    if (length(rmSamples)>0) {
      x = met_data[-rmSamples]
    } else {
      x = met_data
    }
    if (all(is.na(x))) {
      
    } else {
      if (length(x[intersect(which(x>quantile(x, 0.025)), which(x<quantile(x, .975)))]) > 3) {
        x = x[intersect(which(x>quantile(x, 0.025)), which(x<quantile(x, .975)))]
      }
      d = qqnorm(x, plot.it = FALSE);
      x = as.numeric(d$y)
      z = as.numeric(d$x)
      df = data.frame(x=x,z=z)
      t = lm(x~z, data=df)
      mn.est = as.numeric(t$coefficients[1])
      sd.est = as.numeric(t$coefficients[2])
      rm(d,x,z,df,t)
      zscore.data[met,] = (data[met, ]-mn.est)/sd.est
    }
  }
  
  return(zscore.data)
}
