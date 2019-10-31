#' Diffuse Probability P1 from a starting node.
#'
#' Recursively diffuse probability from a starting node based on the connectivity of the background knowledge graph, representing the likelihood that a variable will be
#'         most influenced by a perturbation in the starting node.
#' @param ig_dis - The igraph object associated with the disease+reference trained differential interaction network.
#' @param ig_ref - The igraph object associated with the reference-only trained interaction network.
#' @return ig_pruned - The pruned igraph object of the disease+reference differential interaction network, with reference edges subtracted.
#' @export graph.naivePruning
#' @examples
#' ig_pruned=graph.naivePruning(ig_dis, ig_ref)
graph.naivePruning = function(ig_dis, ig_ref) {
  ig_pruned = ig_dis
  ee = get.edgelist(ig_ref)
  it = 0
  for (e in 1:nrow(ee)) {
    e.id = get.edge.ids(ig_pruned, vp=ee[e,])
    if (e.id != 0) {
      it = it + 1
      ig_pruned = delete.edges(ig_pruned, edges = E(ig_pruned)[[e.id]])
    }
  }
  print(sprintf("%s edges overlapped between reference and disease+reference networks.", it))
  return (ig_pruned)
}






