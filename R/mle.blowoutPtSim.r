#' Module that best explains the patient similarity assigned between a set of patients.
#'
#' @param patientSim - A similarity matrix, where row and columns are patient identifiers.
#' @param data_mx - The matrix that gives the perturbation strength (z-scores) for all variables (columns) for each patient (rows).
#' @param ptIDs - The identifier associated with patient 1's sample.
#' @param ig_pruned - The list of igraph objects associated with the integrated, pruned disease+reference differential interaction networks.
#' @param kmx - The maximum metabolite set size probed when assessing patient similarity.
#' @return ptsim_blowout - An igraph object showing the module blowout describing the similarity between patients in ptIDs.
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
  allmets = unique(unlist(lapply(ig_pruned, function(i) V(i)$name)))
  data_mx = data_mx[which(rownames(data_mx) %in% allmets),]
  # Get all metabolites in top kmx perturbations for all patients in ptIDs.
  # Score node based on how many patient perturbation sets it occurs
  pt_mets = sort(table(as.vector(sapply(ptIDs, function(i) rownames(data_mx)[order(abs(data_mx[,i]), decreasing = TRUE)][1:kmx]))), decreasing = TRUE)
  # Select nodes and edges based on being ranked in top 10% of each scoring.
  selected_nodes = names(pt_mets[which(pt_mets>quantile(pt_mets, 0.50))])
  if (length(grep("x -", selected_nodes)>0)) {
    selected_nodes = selected_nodes[-grep("x -", selected_nodes)]
  }
  
  ptsim_blowout = make_empty_graph(directed=FALSE)
  if (length(selected_nodes)>1) {
    # Rank edges between metabolite nodes by order of strength (absolute value) across all networks in ig_pruned.
    # Keep adding until all nodes are included in a tree
    for (i in 1:length(ig_pruned)) {
      ig = ig_pruned[[i]]
      selected_nodes_ig = selected_nodes[which(selected_nodes %in% V(ig)$name)]
      e_ids = sort(unique(apply(combn(selected_nodes_ig, 2), 2, function(i) get.edge.ids(ig, i))))
      if (any(e_ids==0)) { e_ids= e_ids[-which(e_ids==0)]}
      if (length(e_ids)>0) {
        e_list = get.edgelist(ig)
        e_list2 = e_list[e_ids,]
        if (class(e_list2)=="matrix") {
          rank_edges = e_list2[order(abs(E(ig)$weight[e_ids]), decreasing = TRUE),]
          rank_edgeweights = E(ig)$weight[e_ids][order(abs(E(ig)$weight[e_ids]), decreasing = TRUE)]
          # Rank nodes, display ranked nodes. Node size in blowout is rank of node.
          ind = which(abs(rank_edgeweights)>quantile(abs(E(ig)$weight), 0.95))
          if (length(ind)>1) {
            rank_edges = rank_edges[ind,]
            rank_edges = rank_edges[which(apply(rank_edges, 1, function(i) any(i %in% selected_nodes_ig))),]
          } 
          e_ids2 = sort(apply(rank_edges, 1, function(i) get.edge.ids(ig, i)))
        } else {
          e_ids2 = get.edge.ids(ig, e_list2)
        }
        ig_blowout = subgraph.edges(ig, E(ig)[e_ids2], delete.vertices = TRUE)
        
        vv = V(ig_blowout)$name[which(!(V(ig_blowout)$name %in% V(ptsim_blowout)$name))]
        ptsim_blowout = add.vertices(ptsim_blowout, nv = length(vv), attr=list(name=vv))
        ee = get.edgelist(ig_blowout)
        for (e in 1:nrow(ee)) {
          ptsim_e_id = get.edge.ids(ptsim_blowout, vp=ee[e,])
          if (ptsim_e_id==0) {
            # Edge isn't in ptsim_blowout, add it.
            ptsim_blowout = add.edges(ptsim_blowout, edges = ee[e,], attr=list(weight=E(ig_blowout)$weight[e]))
          } else {
            mx_w = which.max(c(abs(E(ptsim_blowout)$weight[ptsim_e_id]), abs(E(ig_blowout)$weight[e])))
            direction = ifelse(c(E(ptsim_blowout)$weight[ptsim_e_id], E(ig_blowout)$weight[e])[mx_w] <0, -1, 1)
            E(ptsim_blowout)$weight[ptsim_e_id] = direction*max(abs(E(ptsim_blowout)$weight[ptsim_e_id]), abs(E(ig_blowout)$weight[e]))
          }
        }
      }
    }
  } else {
    ptsim_blowout = add.vertices(ptsim_blowout, nv = 1, attr=list(name=selected_nodes))
  }
  #plot.igraph(ptsim_blowout, layout=layout.circle, edge.width=50*abs(E(ptsim_blowout)$weight))
  
  return(ptsim_blowout)
}


