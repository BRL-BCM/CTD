# Loop through all Metabolon Pathways and get a list of unique metabolites.


#' Get All Metabolites In Metabolon's Pathway Knowledgebase
#' @return mets - a character vector of unique metabolites and enzymes found in at least 1 pathway in Metabolon's pathway knowledgebase.
#' @export pathway.ListMetabolites_metabolon
#' @import igraph
#' @examples
#' mets = pathway.ListMetabolites_metabolon()
#' print(mets)
pathway.ListMetabolites_metabolon = function() {
  ig_files = list.files(system.file("extdata/RData/", package="CTD"), pattern = ".RData")
  ig_files = ig_files

  mets = c()
  for (pwy in 1:length(ig_files)) {
    load(system.file(sprintf("extdata/RData/%s", ig_files[pwy]), package="CTD"))
    length(V(ig)$label)
    mets = c(mets, V(ig)$label)
  }
  mets = unique(mets)

  return (mets)
}
