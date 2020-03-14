#' Metabolite set enrichment analysis (MSEA) (using a hypergeometric test) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' @export shiny.getORA_Metabolon
#' @examples
#' pathway.data = shiny.getORA_Metabolon(met.profile, threhold=3, "z-score", NULL)
shiny.getORA_Metabolon = function(input) {
  data = Miller2015[,grep("IEM_", colnames(Miller2015))]
  met.profile = data[,which(colnames(data) %in% input$ptIDs)]
  met.profile = met.profile[which(!(is.na(met.profile)))]
  # TODO: add an input parameter "threshold" to modulate z-score cutoff value.
  perturbed.mets = met.profile[which(abs(met.profile) > 2)]
  nms.perturbed.mets = unique(as.character(unlist(sapply(names(perturbed.mets), function(i) strsplit(i, split=";")))))

  # The size of the population of total possible metabolites to draw from
  population = names(met.profile)
  print(sprintf("Total number of metabolites in profile. = %d.", length(population)))
  paths.hsa = list.dirs(path=system.file("extdata", package="CTD"), full.names = FALSE)
  paths.hsa = paths.hsa[-which(paths.hsa %in% c("", "RData", "allPathways", "MSEA_Datasets", "Pathway_GMTs"))]
  row = 1
  pathway.data = data.frame(Pathway=character(), FDR=numeric(), Pvalue=numeric(), Hits=integer(), Size=integer(), stringsAsFactors = FALSE)
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

  return(pathway.data)
}

