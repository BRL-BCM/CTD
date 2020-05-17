#' shiny.getPathwayIgraph
#'
#' @param input - A list object of parameters (esp. from R shiny app). Required parameters are ptIDs, diagClass and pathwayMapId.
#' @param Pathway.Name - The name of the pathway map for which you want the topological information.
#' @return template.ig - Igraph object of selected pathway map.
#' @export shiny.getPathwayIgraph
#' @usage shiny.getPathwayIgraph(input, Pathway.Name)
#' @import igraph
#' @examples
#' data(Miller2015)
#' # Input is supplied by R shiny app, but you can hard code parameters as a list object, too, to test functionality.
#' input = list()
#' input$ptIDs = colnames(Miller2015)[4]
#' input$diagClass = "paa"
#' input$pathwayMapId = "All"
#' res = shiny.getPathwayIgraph(input, Miller2015)
#' # Returns a blank template for selected pathway, node labels and node types.

shiny.getPathwayIgraph = function(input, Pathway.Name) {
  if (is.null( Pathway.Name)) Pathway.Name = gsub(" ", "-", input$pathwayMapId)

  if (Pathway.Name=="All") {
    load(system.file("extdata/RData/allPathways.RData", package="CTD"))
    V(ig)$label[which(V(ig)$label %in% c("DSGEGDFXAEGGGVR", "Dsgegdfxaegggvr"))] = ""
    Pathway.Name = "allPathways"
  } else {
    load(system.file(sprintf("extdata/RData/%s.RData", Pathway.Name), package="CTD"))
  }
  template.ig = ig

  # Load id to node types mappings
  nodeType = read.table(system.file(sprintf("extdata/%s/Type-%s.txt", Pathway.Name, Pathway.Name), package="CTD"), header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = suppressWarnings(as.numeric(tmp.nms))
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeType = as.character(tmp)
  names(nodeType) = tmp.nms
  nodeType = nodeType[which(names(nodeType) %in% V(template.ig)$name)]

  # node.labels = vector("character", length = length(V(template.ig)$name))
  node.labels=V(template.ig)$label
  node.types = vector("character", length = length(V(template.ig)$name))
  for (n in 1:length(V(template.ig)$name)) {
    # node.labels[n] = URLdecode(as.character(nodeDisplayNames[V(template.ig)$name[n]]))
    node.types[n] = as.character(nodeType[V(template.ig)$name[n]])
  }

  return(list(template.ig=template.ig, nodeDisplayNames=node.labels, nodeType=node.types))
}
