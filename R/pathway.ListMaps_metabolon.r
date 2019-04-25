# Loop through all Metabolon Pathways and get a list of unique metabolites.


#' Get All Metabolites In Metabolon's Pathway Knowledgebase
#' @return Pathway.Knowledgebase.Path - local path to extdata folder containing metabolon's pathway knowledgebase files.
#' @export pathway.ListMaps_metabolon
#' @import igraph
#' @examples
#' metabolon_knowledgebase_path = sprintf("%s/extdata", find.package("CTD"))
#' pwys = pathway.ListMaps_metabolon(metabolon_knowledgebase_path)
#' print(pwys)
pathway.ListMaps_metabolon = function(Pathway.Knowledgebase.Path) {
  ig_files = list.files(sprintf("%s/extdata/RData/", Pathway.Knowledgebase.Path), pattern = ".RData")
  ig_files = ig_files
  pwys = unlist(sapply(ig_files, function(i) unlist(strsplit(i, split=".RData"))[1]))
  pwys = as.character(pwys)
  pwys = gsub("-", " ", pwys)
  pwys[which(pwys=="allPathways")] = "All Pathways"

  return (pwys)
}
