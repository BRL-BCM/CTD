#' Visualize the confusion matrix using nearest neighbor as classification model.
#'
#' @param simMat - The patient similarity matrix.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of simMat. 
#'                    The names of this vector are patient IDs, and values are diagnostic labels.
#' @param k - The number of dimension you want to plot your data using multi-dimensional scaling.
#' @param diag - The diagnosis associated with positive controls in your data.
#' @return p - a plotly scatter plot colored by provided diagnostic labels.
#' @export plot.knnSim
#' @examples
#' # Look at main_CTD.r script for full analysis script: https://github.com/BRL-BCM/CTD.
#' # if you have diagnostic labels associated with the colnames(data_mx), send them using diagnoses parameter
#' p = plot.knnSim(patientSim, diagnoses, k=2, diag="diseased")
#' p
plot.knnSim = function(patientSim, diagnoses, k, diag) {
  # Add a GREEN edge between patient nodes if k nearest neighbor is correct diagnosis (either TP or TN)
  # Add a RED edge between patient nodes if k nearest neighbor is incorrect diagnosis (either FP or FN)
  ig = make_empty_graph(n=ncol(patientSim), directed=TRUE)
  V(ig)$name = colnames(patientSim)
  for (pt1 in 1:length(diagnoses)) {
    diag_pt1 = diagnoses[pt1]
    ind = sort(patientSim[pt1,-pt1], decreasing = FALSE)
    diag_pt_ind = diagnoses[which(names(diagnoses)==names(ind[1]))]
    if (diag_pt_ind==diag && diag_pt1==diag) {
      # True positive
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="green", lty=1))
    } else if (diag_pt_ind!=diag && diag_pt1!=diag) {
      # True negative
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="green", lty=1))
    } else {
      # Mistake! Either false positive or false negative
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="red", lty=1))
    }
  }
  V(ig)$label = rep("", length(V(ig)$name))
  tt = table(diagnoses)
  grps = list()
  grps[1] = names(diagnoses)[which(diagnoses==diag)]
  names(grps[1]) = diag
  grps[2] = names(diagnoses)[which(diagnoses!=diag)]
  names(grps[2]) = "controls"
  p = plot(ig, mark.groups = grps, mark.col = c("grey", "black"),
       layout=layout.circle, edge.width=3, edge.arrow.size=0.5) + 
    legend("topright", legend=c(diag, "control"), fill=c("grey", "black"))
  
  return(p)
}






