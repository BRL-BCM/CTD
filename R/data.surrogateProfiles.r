#' Surrogate profiles
#'
#' Fill in a data matrix with low n, high p with surrogate profiles.
#' @param data - Data matrix with observations as rows, features as columns.
#' @param sd - The level of variability (standard deviation) around each feature's
#'             mean you want to add in surrogate profiles.
#' @return data_mx - Data matrix with added surrogate profiles.
#' @export data.surrogateProfiles
data.surrogateProfiles = function(data, sd=1, useMnProfile=FALSE, ref_data=NULL) {

  if (useMnProfile==TRUE) {
    mn_data = apply(data, 1, function(i) mean(na.omit(i)))
    names(mn_data) = rownames(data)
    data = as.matrix(mn_data)
    colnames(data) = "mn_data"
    # Generate disease surrogates
    rpt = ceiling(nrow(data)/ncol(data)/2)
    data_surr = matrix(NA, nrow=nrow(data), ncol=ncol(data)+ncol(data)*(rpt+1))
    c_col = ncol(data)+1
    for (pt in 1:ncol(data)) {
      data_surr[,pt] = data[,pt]
      for (rrpt in 1:(rpt+1)) {
        rr = rnorm(nrow(data), mean=0, sd=sd)
        data_surr[,c_col] = as.numeric(data[,pt])+rr
        c_col = c_col+1
      }
    }
    colnames(data_surr) = c(colnames(data), sprintf("disease_surr%d", 1:(ncol(data_surr)-ncol(data))))
    rownames(data_surr) = rownames(data)
    dim(data_surr)
  } else {
    # Generate disease surrogates
    rpt = ceiling(nrow(data)/ncol(data)/2)
    data_surr = matrix(NA, nrow=nrow(data), ncol=ncol(data)+ncol(data)*(rpt+1))
    c_col = ncol(data)+1
    for (pt in 1:ncol(data)) {
      data_surr[,pt] = data[,pt]
      for (rrpt in 1:(rpt+1)) {
        rr = rnorm(nrow(data), mean=0, sd=sd)
        data_surr[,c_col] = as.numeric(data[,pt])+rr
        c_col = c_col+1
      }
    }
    colnames(data_surr) = c(colnames(data), sprintf("disease_surr%d", 1:(ncol(data_surr)-ncol(data))))
    rownames(data_surr) = rownames(data)
    dim(data_surr)
  }

  if (!is.null(ref_data)) {
    ref_data = ref_data[which(rownames(ref_data) %in% rownames(data)), ]
    # Generate control surrogates
    numSurr = ceiling(nrow(data)/2)
    if (numSurr> ncol(ref_data)) {
      numSurr = numSurr - ncol(ref_data)
      control_surr = matrix(0, nrow=nrow(data), ncol=numSurr)
      for (pt in 1:numSurr) {
        control_surr[,pt] = rnorm(nrow(data), mean=0, sd=sd)
      }
      control_surr = cbind(ref_data, control_surr)
    } else {
      control_surr = ref_data[,sample(1:ncol(ref_data), numSurr)]
    }
    colnames(control_surr) = sprintf("control_surr%d", 1:numSurr)
    rownames(control_surr) = rownames(ref_data)
    dim(control_surr)
  } else {
    # Generate control surrogates
    numSurr = ceiling(nrow(data)/2)
    control_surr = matrix(0, nrow=nrow(data), ncol=numSurr)
    for (pt in 1:numSurr) {
      control_surr[,pt] = rnorm(nrow(data), mean=0, sd=sd)
    }
    colnames(control_surr) = sprintf("control_surr%d", 1:numSurr)
    rownames(control_surr) = rownames(data)
    dim(control_surr)
  }


  data_mx = cbind(data_surr, control_surr)
  rmThese = c()
  for (r in 1:nrow(data_mx)) {
    if (all(is.na(as.numeric(data_mx[rownames(data_mx)[r],])))) {
      rmThese = c(rmThese, r)
    } else {
      data_mx[r, which(is.na(data_mx[r,]))] = min(na.omit(as.numeric(control_surr[rownames(data_mx)[r],])))
    }
  }
  if (length(rmThese)>0) {
    data_mx = data_mx[-rmThese,]
  }
  any(is.na(data_mx)) # Should be FALSE
  any(is.infinite(unlist(data_mx))) # Should be FALSE

  # remove any metabolite with no variation (sd=0)
  var.met = apply(data_mx, 1, sd)
  if (length(which(var.met==0))>0) {
    data_mx = data_mx[-which(var.met==0),]
  }
  rownames(data_mx) = tolower(rownames(data_mx))
  data_mx = apply(data_mx, c(1,2), as.numeric)
  dim(data_mx)

  return(data_mx)
}
