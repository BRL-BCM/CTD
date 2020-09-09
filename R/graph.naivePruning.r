#' Network pruning for disease-specific network determination
#' 
#' Prune edges from a disease+control "differential" network that also 
#' occur in the control-only network.
#' @param ig_dis - The igraph object associated with the
#'                 disease+reference trained differential network.
#' @param ig_ref - The igraph object associated with the reference-only
#'                 trained interaction network.
#' @return ig_pruned - The pruned igraph object of the disease+reference
#' differential network, with reference edges subtracted.
#' @import igraph
#' @export graph.naivePruning
#' @examples
#' # Generate a 100 node "disease-control" network
#' adj_mat=matrix(0, nrow=100, ncol=100)
#' rows = sample(seq_len(100), 50, replace=TRUE)
#' cols = sample(seq_len(100), 50, replace=TRUE)
#' for (i in rows) {for (j in cols){adj_mat[i,j]=rnorm(1,0,1)}}
#' colnames(adj_mat)=sprintf("Metabolite%d", seq_len(100))
#' ig_dis = graph.adjacency(adj_mat, mode="undirected", weighted=TRUE)
#' # Generate a 100 node reference "control-only" network
#' adj_mat2=matrix(0, nrow=100, ncol=100)
#' rows2 = sample(seq_len(100), 50, replace=TRUE)
#' cols2 = sample(seq_len(100), 50, replace=TRUE)
#' for (i in rows2) {for (j in cols2){adj_mat2[i,j]=rnorm(1,0,1)}}
#' colnames(adj_mat2)=sprintf("Metabolite%d", seq_len(100))
#' ig_ref = graph.adjacency(adj_mat2, mode="undirected", weighted=TRUE)
#' ig_pruned=graph.naivePruning(ig_dis, ig_ref)
graph.naivePruning = function(ig_dis, ig_ref) {
    ee=get.edgelist(ig_ref)
    ee=ee[which(apply(ee,1,function(i) all(i %in% V(ig_dis)$name))),]
    if (is.null(nrow(ee))) { ee = t(ee) }
    it=0
    for (e in seq_len(nrow(ee))) {
        e.id.dis=get.edge.ids(ig_dis, vp=ee[e,])
        e.id.ref=get.edge.ids(ig_ref, vp=ee[e,])
        if (e.id.dis!=0) {
            isSameDirection=ifelse((((E(ig_dis)$weight[e.id.dis]<0) &&
                                        (E(ig_ref)$weight[e.id.ref]<0)) || 
                                    ((E(ig_dis)$weight[e.id.dis]>0) &&
                                        (E(ig_ref)$weight[e.id.ref]>0))),
                                    TRUE,FALSE)
            if (isSameDirection && (abs(E(ig_ref)$weight[e.id.ref])<
                                    abs(E(ig_dis)$weight[e.id.dis]))) {
                new_w=E(ig_dis)$weight[e.id.dis]-E(ig_ref)$weight[e.id.ref]
                E(ig_dis)$weight[e.id.dis]=new_w
                it=it + 1
            }
            if (isSameDirection && (abs(E(ig_ref)$weight[e.id.ref]) >=
                                    abs(E(ig_dis)$weight[e.id.dis]))) {
                ig_dis=delete.edges(ig_dis, edges = E(ig_dis)[[e.id.dis]])
                it=it + 1
            }
        }
    }
    print(sprintf("%s edges were modified in the dis+ref net.", it))
    return (ig_dis)
}