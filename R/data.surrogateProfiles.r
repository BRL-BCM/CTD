#' Surrogate profiles
#'
#' Fill in a data matrix rank, when your data is low n, high p. Fill in rank with surrogate profiles.
#' @param data - Data matrix with observations as rows, features as columns.
#' @param sd - The level of variability (standard deviation) around each feature's
#'             mean you want to add in surrogate profiles.
#' @param useMnDiseaseProfile - Boolean. For disease cohorts not showing homogeneity, mean
#'             across disease profiles and generate disease surrogates around this mean.
#' @param addHealthyControls - Boolean. Add healthy control profiles to data?
#' @return data_mx - Data matrix with added surrogate profiles.
#' @export data.surrogateProfiles
data.surrogateProfiles = function(data, sd=1, useMnDiseaseProfile=FALSE, addHealthyControls=TRUE, ref_data=NULL) {
  ref_data = ref_data[which(rownames(ref_data) %in% rownames(data)),]
  data = data[which(rownames(data) %in% rownames(ref_data)), ]

  # Generate disease surrogates first
  if (useMnDiseaseProfile==TRUE) {
    mn_data = apply(data, 1, function(i) mean(na.omit(i)))
    names(mn_data) = rownames(data)
    data = as.matrix(mn_data)
    colnames(data) = "mn_data"
    # Generate disease surrogates
    if (addHealthyControls) {
      rpt = ceiling(nrow(data)/ncol(data)/2)
    } else {
      rpt = ceiling(nrow(data)/ncol(data))
    }
    data_surr = matrix(NA, nrow=nrow(data), ncol=ncol(data)+ncol(data)*rpt)
    c_col = ncol(data)+1
    for (pt in 1:ncol(data)) {
      data_surr[,pt] = data[,pt]
      for (rrpt in 1:rpt) {
        rr = rnorm(nrow(data), mean=0, sd=sd)
        data_surr[,c_col] = as.numeric(data[,pt])+rr
        c_col = c_col+1
      }
    }
    colnames(data_surr) = c(colnames(data), sprintf("disease_surr%d", 1:(ncol(data_surr)-ncol(data))))
    rownames(data_surr) = rownames(data)
    dim(data_surr)
  } else {
    # Generate disease surrogates for each "flavor" of disease profiles.
    if (addHealthyControls) {
      rpt = ceiling(nrow(data)/ncol(data)/2)
    } else {
      rpt = ceiling(nrow(data)/ncol(data))
    }
    data_surr = matrix(NA, nrow=nrow(data), ncol=ncol(data)+ncol(data)*rpt)
    c_col = ncol(data)+1
    for (pt in 1:ncol(data)) {
      data_surr[,pt] = data[,pt]
      for (rrpt in 1:rpt) {
        rr = rnorm(nrow(data), mean=0, sd=sd)
        data_surr[,c_col] = as.numeric(data[,pt])+rr
        c_col = c_col+1
      }
    }
    colnames(data_surr) = c(colnames(data), sprintf("disease_surr%d", 1:(ncol(data_surr)-ncol(data))))
    rownames(data_surr) = rownames(data)
    dim(data_surr)
  }

  # Next, if addHealthyControls is TRUE, generate reference profile surrogates.
  if (addHealthyControls) {
    numSurr = ceiling(nrow(data)/2)
    if (numSurr > ncol(ref_data)) {
      # Generate control surrogates from real control samples.
      rpt = ceiling(nrow(data)/ncol(ref_data)/2)
      control_surr = matrix(NA, nrow=nrow(ref_data), ncol=ncol(ref_data)+ncol(ref_data)*rpt)
      c_col = ncol(ref_data)+1
      for (pt in 1:ncol(ref_data)) {
        control_surr[,pt] = ref_data[,pt]
        for (rrpt in 1:rpt) {
          rr = rnorm(nrow(ref_data), mean=0, sd=sd)
          control_surr[,c_col] = as.numeric(ref_data[,pt])+rr
          c_col = c_col+1
        }
      }
      colnames(control_surr) = c(colnames(ref_data), sprintf("control_surr%d", 1:(ncol(control_surr)-ncol(ref_data))))
      rownames(control_surr) = rownames(ref_data)
      dim(control_surr)
    } else {
      ind = sample(1:ncol(ref_data), numSurr)
      control_surr = ref_data[,ind]
      colnames(control_surr) = colnames(control_surr)[ind]
      rownames(control_surr) = rownames(ref_data)
      dim(control_surr)
    }
    data_mx = cbind(data_surr, control_surr)
  } else {
    data_mx = data_surr
  }

  # Remove metabolites that are NA for all samples in data_mx
  rmThese = c()
  for (r in 1:nrow(data_mx)) {
    if (all(is.na(as.numeric(data_mx[rownames(data_mx)[r],])))) {
      rmThese = c(rmThese, r)
    } else {
      data_mx[r, which(is.na(data_mx[r, ]))] = min(na.omit(as.numeric(ref_data[rownames(data_mx)[r],])))
    }
  }
  if (length(rmThese) > 0) { data_mx = data_mx[-rmThese, ] }
  any(is.na(data_mx))
  any(is.infinite(unlist(data_mx)))

  var.met = apply(data_mx, 1, sd)
  if (length(which(var.met == 0)) > 0) { data_mx = data_mx[-which(var.met == 0), ] }
  rownames(data_mx) = tolower(rownames(data_mx))
  data_mx = apply(data_mx, c(1, 2), as.numeric)
  dim(data_mx)

  return(data_mx)
}
