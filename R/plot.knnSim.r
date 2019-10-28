#' Visualize the confusion matrix using nearest neighbor as classification model.
#'
#' @param patientSim - The patient similarity matrix.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of patientSim.
#'                    The names of this vector are patient IDs, and values are diagnostic labels.
#' @param diag - The diagnosis associated with positive controls in your data.
#' @return p - a plotly scatter plot colored by provided diagnostic labels.
#' @export plot.knnSim
#' @usage plot.knnSim(patientSim, diagnoses, diag)
#' @examples
#' # if you have diagnostic labels associated with the colnames(data_mx), send them using diagnoses parameter
#' p = plot.knnSim(patientSim, diagnoses, diag="diseased")
#' p
plot.knnSim = function(patientSim, diagnoses, k, diag) {
  diagnoses = diagnoses[colnames(patientSim)]
  # Add a GREEN edge between patient nodes if k nearest neighbor is correct diagnosis (either TP or TN)
  # Add a RED edge between patient nodes if k nearest neighbor is incorrect diagnosis (either FP or FN)
  tp = 0
  fp = 0
  tn = 0
  fn = 0
  ig = make_empty_graph(n=ncol(patientSim), directed=TRUE)
  V(ig)$name = colnames(patientSim)
  for (pt1 in 1:length(diagnoses)) {
    diag_pt1 = diagnoses[pt1]
    ind = sort(patientSim[pt1,-pt1], decreasing = FALSE)
    ind = ind[which(ind==min(ind))]
    diag_pt_ind = diagnoses[which(names(diagnoses) %in% names(ind))]
    if (any(diag_pt_ind==diag) && diag_pt1==diag) {
      # True positive
      tp = tp + 1
      ind = ind[which(diag_pt_ind==diag)]
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="green", lty=1))
    } else if (diag_pt_ind!=diag && diag_pt1!=diag) {
      # True negative
      tn = tn + 1
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="green", lty=1))
    } else if (diag_pt_ind==diag && diag_pt1!=diag) {
      # Mistake!
      fp = fp + 1
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="red", lty=1))
    } else {
      fn = fn + 1
      ig = add.edges(ig, edges=c(colnames(patientSim)[pt1], names(ind[1])), attr=list(color="red", lty=1))
    }
  }
  print(sprintf("Tp = %d, Tn= %d, Fp = %d, Fn=%d", tp, tn, fp, fn))
  sens = tp / (tp+fn)
  spec = tn / (tn+fp)
  print(sprintf("Sens = %.2f, Spec= %.2f", sens, spec))

  V(ig)$label = rep("", length(V(ig)$name))
  grps = list()
  grps[[1]] = names(diagnoses)[which(diagnoses==diag)]
  grps[[2]] = names(diagnoses)[which(diagnoses!=diag)]
  names(grps) = names(table(diagnoses))
  plot(ig, mark.groups = grps, mark.col = c("white", "black"), layout=layout.circle, edge.width=3, edge.arrow.size=0.5) #+
    #legend("topright", legend=names(grps), fill=c("white", "black"))

  #tt = table(diagnoses)
  #grps = list()
  #for (i in 1:length(tt)) {
  #  grps[[i]] = names(diagnoses)[which(diagnoses==names(tt)[i])]
  #}
  #names(grps) = names(tt)
  #p = plot(ig, mark.groups = grps, mark.col = c("white", "light grey", "dark grey", "black")[1:length(tt)], layout=layout.circle, edge.width=3, edge.arrow.size=0.5) +
  #  legend("topright", legend=names(grps), fill=c("white", "light grey", "dark grey", "black")[1:length(tt)])

  return(0)
}






