#' View patient clusters using multi-dimensional scaling.
#'
#' This function plots the provided patient similarity matrix in a lower dimensional space
#' using multi-dimensional scaling, which is well suited for similarity metrics.
#' @param simMat - The patient similarity matrix.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of simMat.
#' @param k - The number of dimension you want to plot your data using multi-dimensional scaling.
#' @export plot.mdsSim
#' @examples
#' plot.mdsSim(simMat, path)
plot.mdsSim = function(simMat, diagnoses, k) {
  if (!(k %in% c(2,3))) {
    print("K must be either 2-dimensions or 3-dimensions.")
    return(0)
  }

  if (is.null(diagnoses)) {
    print("To view patient clusters, please provide clinical labels, corresponding to each patient in the provided similarity matrix.")
    return(0)
  }

  fitSim = cmdscale(patientSim, eig=FALSE, k=k)
  x = round(fitSim[,1], 2)
  y = round(fitSim[,2], 2)
  if (k==3) {
    z = round(fitSim[,3], 2)
    df = data.frame(x=x, y=y, z=z, color=diagnoses, label=colnames(patientSim))
    p = plot_ly(df, x=~x, y=~y, z=~z, color=~color, text=~label, marker = list(size = 20))
  } else {
    df = data.frame(x=x, y=y, color=diagnoses, label=colnames(patientSim))
    p = plot_ly(df, x=~x, y=~y, color=~color, text=~label, marker = list(size = 20))
  }

  return(p)
}


