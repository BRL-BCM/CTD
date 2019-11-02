# Loop through all Metabolon Pathways and get a list of unique metabolites.


#' Get All Metabolites In Metabolon's Pathway Knowledgebase
#' @return pwys - List of pathway maps curated by Metabolon's Metabolync.
#' @export pathway.ListMaps_metabolon
#' @import igraph
#' @examples
#' pwys = pathway.ListMaps_metabolon()
#' print(pwys)
pathway.ListMaps_metabolon = function() {
  ig_files = list.files(system.file("extdata/RData/", package="CTD"), pattern = ".RData")
  ig_files = ig_files
  pwys = unlist(sapply(ig_files, function(i) unlist(strsplit(i, split=".RData"))[1]))
  pwys = as.character(pwys)
  pwys = gsub("-", " ", pwys)
  pwys[which(pwys=="allPathways")] = "All Pathways"

  return (pwys)
}
