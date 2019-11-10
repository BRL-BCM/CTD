#' Patient similarity using mutual information MLE metric of patients' most modular, perturbed subsets.
#'
#' This function calculates the universal distance between patients, using a mutual information metric, where self-information comes from the minimum encoding length of each patient's encoded modular perturbations in the background knowledge graph.
#' @param patientSim - A similarity matrix, where row and columns are patient identifiers.
#' @param data_mx - The matrix that gives the perturbation strength (z-scores) for all variables (columns) for each patient (rows).
#' @param ptIDs - The identifier associated with patient 1's sample.
#' @param ig_pruned - The igraph object associated with the pruned disease+reference differential interaction network.
#' @param kmx - The maximum metabolite set size probed when assessing patient similarity.
#' @return p - A ggplot2 object showing the module blowout describing the similarity between patients in ptIDs.
#' @export mle.blowoutSim
#' @examples
#' require(CTD)
#' data(Thistlethwaite2019)
#' data_mx = as.matrix(data_mx)
#' data_mx = suppressWarnings(apply(data_mx, c(1,2), as.numeric))
#' data_mx = data_mx[,-c(1,2,3,4,5,6,7,8)]
#' # Load your background network, ig_pruned and your computed patientSim matrix
#' kmns.clust = kmeans(patientSim, centers=4)
#' table(kmns.clust$cluster)
#' ptIDs = names(kmns.clust$cluster[which(kmns.clust$cluster==1)])
#' ptsim_blowout = mle.blowoutSim(patientSim, data_mx, ptIDs, ig_pruned, kmx=15)
#' plot.igraph(ptsim_blowout, layout=layout.circle, edge.width=50*abs(E(ptsim_blowout)$weight))
mle.blowoutSim = function(patientSim, data_mx, ptIDs, ig_pruned, kmx) {
  data_mx = data_mx[which(rownames(data_mx) %in% V(ig_pruned)$name),]
  # Get all metabolites in top kmx perturbations for all patients in ptIDs.
  # Score node based on how many patient perturbation sets it occurs
  pt_mets = sort(table(as.vector(sapply(ptIDs, function(i) rownames(data_mx)[order(abs(data_mx[,i]), decreasing = TRUE)][1:kmx]))), decreasing = TRUE)
  # Select nodes and edges based on being ranked in top 10% of each scoring.
  selected_nodes = names(pt_mets[which(pt_mets>quantile(pt_mets, 0.50))])
  
  # Rank edges between metabolite nodes by order of strength (absolute value)
  # Keep adding until all nodes are included in a tree
  e_ids = sort(unique(apply(combn(names(pt_mets), 2), 2, function(i) get.edge.ids(ig_pruned, i))))
  if (any(e_ids==0)) { e_ids= e_ids[-which(e_ids==0)]}
  e_list = get.edgelist(ig_pruned)
  e_list2 = e_list[e_ids,]
  rank_edges = e_list2[order(abs(E(ig_pruned)$weight[e_ids]), decreasing = TRUE),]
  rank_edgeweights = E(ig_pruned)$weight[e_ids][order(abs(E(ig_pruned)$weight[e_ids]), decreasing = TRUE)]
  # Rank nodes, display ranked nodes. Node size in blowout is rank of node.
  rank_edges = rank_edges[which(abs(rank_edgeweights)>quantile(abs(rank_edgeweights), 0.5)),]
  rank_edges = rank_edges[which(apply(rank_edges, 1, function(i) any(i %in% selected_nodes))),]
  
  e_ids2 = sort(apply(rank_edges, 1, function(i) get.edge.ids(ig_pruned, i)))
  
  ptsim_blowout = subgraph.edges(ig_pruned, E(ig_pruned)[e_ids2], delete.vertices = TRUE)
  #plot.igraph(ptsim_blowout, layout=layout.circle, edge.width=50*abs(E(ptsim_blowout)$weight))
  
  return(ptsim_blowout)
}