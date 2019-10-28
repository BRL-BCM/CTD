#' View patient clusters using multi-dimensional scaling.
#'
#' This function plots the provided patient similarity matrix in a lower dimensional space
#' using multi-dimensional scaling, which is well suited for similarity metrics.
#' @param patientSim - The patient similarity matrix.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of patientSim.
#' @param k - The number of dimension you want to plot your data using multi-dimensional scaling.
#' @param diag - The diagnosis associated with positive controls in your data.
#' @return p - a plotly scatter plot colored by provided diagnostic labels.
#' @export plot.mdsSim
#' @usage plot.mdsSim(patientSim, diagnoses, k, diag)
#' @examples
#' # if you have diagnostic labels associated with the colnames(data_mx), send them using diagnoses parameter
#' p = plot.mdsSim(patientSim, diagnoses, k=2, diag="diseased")
#' p
#' p = plot.mdsSim(patientSim, diagnoses, k=3, diag="diseased")
#' p
plot.mdsSim = function(patientSim, diagnoses, k, diag) {
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

  # Do K-means here
  #dis.centroid = apply(patientSim[which(diagnoses==diag),], 2, mean)
  #negCntl.centroid = apply(patientSim[which(diagnoses=="negCntl"),], 2, mean)
  #kmns = kmeans(patientSim, centers = rbind(dis.centroid, negCntl.centroid))
  #dcntd = diagnoses[which(colnames(patientSim) %in% names(kmns$cluster[which(kmns$cluster==1)]))]
  #ncntd = diagnoses[which(colnames(patientSim) %in% names(kmns$cluster[which(kmns$cluster==2)]))]

  #rownames(df) = df[,"label"]
  #min(df[dcntd,"x"])
  #max(df[dcntd,"x"])
  #min(df[dcntd,"y"])
  #max(df[dcntd,"y"])

  #min(df[ncntd,"x"])
  #max(df[ncntd,"x"])
  #min(df[ncntd,"y"])
  #max(df[ncntd,"y"])

  # https://plot.ly/r/shapes/
  # Circles are centered around ((x0+x1)/2, (y0+y1)/2))
  #p = layout(p, title = sprintf("Kmeans Clustering of Patient Distances from MolPhenoMatch"),
  #           shapes = list(
  #             list(type = 'circle',
  #                  xref = 'x', x0 = .2, x1 = .7,
  #                  yref = 'y', y0 = 20, y1 = 3000,
  #                  fillcolor = 'rgb(50, 20, 90)', line = list(color = 'rgb(50, 20, 90)'),
  #                  opacity = 0.2),
  #             list(type = 'circle',
  #                  xref = 'x', x0 = .75, x1 = 1.5,
  #                  yref = 'y', y0 = 2500, y1 = 7500,
  #                  fillcolor = 'rgb(30, 100, 120)', line = list(color = 'rgb(30, 100, 120)'),
  #                  opacity = 0.2)))

  return(p)
}
