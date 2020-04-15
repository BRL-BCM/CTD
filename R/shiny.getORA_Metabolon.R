#' Metabolite set enrichment analysis (MSEA) (using a hypergeometric test) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' @export shiny.getORA_Metabolon
#' @examples
#' pathway.data = shiny.getORA_Metabolon(met.profile, threhold=3, "z-score", NULL)
shiny.getORA_Metabolon = function(input) {
  data = .GlobalEnv$data_zscore
  tmp = rownames(data)
  met.profile = apply(as.matrix(data[,which(colnames(data) %in% input$ptIDs)]), 1, mean)
  names(met.profile) = tmp
  met.profile = met.profile[which(!(is.na(met.profile)))]
  # TODO: add an input parameter "threshold" to modulate z-score cutoff value.
  perturbed.mets = met.profile[which(abs(met.profile) > 2)]
  nms.perturbed.mets = unique(as.character(unlist(sapply(names(perturbed.mets), function(i) strsplit(i, split=";")))))

  # The size of the population of total possible metabolites to draw from
  population = names(met.profile)
  print(sprintf("Total number of metabolites in profile. = %d.", length(population)))
  paths.hsa = list.dirs(path=system.file("extdata", package="CTD"), full.names = FALSE)
  paths.hsa = paths.hsa[-which(paths.hsa %in% c("", "RData", "allPathways", "MSEA_Datasets"))]
  row = 1
  pathway.data = data.frame(Pathway=character(), Size=integer(), Hits=integer(), FDR=numeric(), Pvalue=numeric(), stringsAsFactors = FALSE)
  for (pathway in 1:length(paths.hsa)) {
    load(system.file(sprintf("extdata/RData/%s.RData", paths.hsa[pathway]), package="CTD"))

    pathway.compounds = tolower(V(ig)$label[which(V(ig)$shape=="circle")])
    pathCompIDs = unique(tolower(pathway.compounds[which(pathway.compounds %in% population)]))
    print(sprintf("Total number of metabolites in profile also in current Metabolon pathway. = %d.", length(pathCompIDs)))

    # q (sample successes), m (population successes), n (population failures), k (sample size)
    sampleSuccesses = length(which(nms.perturbed.mets %in% pathCompIDs))
    populationSuccesses = length(intersect(pathCompIDs, population))
    N = length(population)
    populationFailures=N-populationSuccesses
    numDraws=length(perturbed.mets)
    if (populationSuccesses>0) {
      pathway.data[row, "Pathway"] = paths.hsa[pathway]
      pathway.data[row, "Pvalue"] = phyper(q=sampleSuccesses-1, m=populationSuccesses, n=populationFailures, k=numDraws, lower.tail=FALSE)
      pathway.data[row, "Hits"] = sampleSuccesses
      pathway.data[row, "Size"] = populationSuccesses
      row = row + 1
    }
  }
  pathway.data[,"FDR"] = p.adjust(pathway.data[,"Pvalue"], method="fdr")
  # Sort by FDR, then Pvalue
  pathway.data = pathway.data[order(pathway.data[,"FDR"], pathway.data[,"Pvalue"]),]

  # Then, only return rows that have Pvalue < 0.25
  pathway.data = pathway.data[which(pathway.data[,"Pvalue"]<0.25),]
  
  pathway.data$FDR = signif(pathway.data$FDR, digits=3)
  pathway.data$Pvalue = signif(pathway.data$Pvalue, digits=3)
  
  return(pathway.data)
}

