#' Generate heatmap plot of patient similarity matrix.
#'
#' This function plots a heatmap of a patient similarity matrix.
#' @param patientSim - The patient similarity matrix.
#' @param path - The filepath to a directory in which you want to store the .png file.
#' @param diagnoses - A character vector of diagnostic labels associated with the rownames of patientSim.
#' @export plot.hmSim
#' @usage plot.hmSim(patientSim, path=getwd(), diagnoses=NULL)
#' @examples
#' plot.hmSim(patientSim, path=getwd(), diagnoses)
plot.hmSim = function(patientSim, path, diagnoses=NULL) {
  if (is.null(diagnoses)) {
    png(sprintf("%s/ptSimilarity.png", path), 3000, 1000)
    heatmap.2(x=patientSim,
              dendrogram="both",
              Rowv=TRUE,
              Colv=TRUE,
              cexRow=0.75,cexCol=1, margins=c(12,12),
              trace="none", key=TRUE,
              col=bluered,
              notecol="black")
    dev.off()
  } else {
    d = diagnoses
    png(sprintf("%s/ptSimilarity.png", path), 3000, 1000)
    heatmap.2(x=patientSim,
              dendrogram="both",
              Rowv=TRUE,
              Colv=TRUE,
              ColSideColors = c(brewer.pal(12, "Set3"), brewer.pal(9, "BrBG"))[as.numeric(as.factor(as.character(d)))],
              RowSideColors = c(brewer.pal(12, "Set3"), brewer.pal(9, "BrBG"))[as.numeric(as.factor(as.character(d)))],
              cexRow=0.75,cexCol=1, margins=c(12,12),
              trace="none", key=TRUE,
              col=bluered,
              notecol="black")
    legend("left", legend=unique(sort(as.character(d))), fill=c(brewer.pal(12, "Set3"), brewer.pal(9, "BrBG")), cex=2)
    dev.off()
  }
  return(0)
}


