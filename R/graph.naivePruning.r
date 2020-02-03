#' Pruning edges from disease differential network that also occur in reference-only network.
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
  ee = get.edgelist(ig_ref)
  ee = ee[which(apply(ee, 1, function(i) all(i %in% V(ig_dis)$name))),]
  it = 0
  for (e in 1:nrow(ee)) {
    e.id.dis = get.edge.ids(ig_dis, vp=ee[e,])
    e.id.ref = get.edge.ids(ig_ref, vp=ee[e,])
    if (e.id.dis != 0) {
      isSameDirection = ifelse((((E(ig_dis)$weight[e.id.dis]<0) && (E(ig_ref)$weight[e.id.ref]<0)) || 
                               ((E(ig_dis)$weight[e.id.dis]>0) && (E(ig_ref)$weight[e.id.ref]>0))), TRUE, FALSE)
      
      if (isSameDirection && (abs(E(ig_ref)$weight[e.id.ref]) < abs(E(ig_dis)$weight[e.id.dis]))) {
        E(ig_dis)$weight[e.id.dis] = E(ig_dis)$weight[e.id.dis] - E(ig_ref)$weight[e.id.ref]
        it = it + 1
      }
      if (isSameDirection && (abs(E(ig_ref)$weight[e.id.ref]) > abs(E(ig_dis)$weight[e.id.dis]))) {
        ig_dis = delete.edges(ig_dis, edges = E(ig_dis)[[e.id.dis]])
        it = it + 1
      }
    }
  }
  print(sprintf("%s edges were modified in the disease+reference network.", it))
  return (ig_dis)
}






