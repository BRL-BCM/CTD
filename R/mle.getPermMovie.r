#' Capture the movement of the adaptive walk of the diffusion probability method.
#'
#' Make a movie of the adaptive walk the diffusion probability method makes in search of a given patient's perturbed variables.
#' @param patient.sig.nodes - The subset of variables, S, in a background graph, G. These are the perturbed molecular phenotype variables in a given patient's profile.
#' @param patient - The column number in data.pvals associated with the patient being examined.
#' @param ig - The igraph object associated with the background knowledge graph.
#' @param movie - If you want to make a movie, set to TRUE. This will produce a set of still images that you can stream together to make a movie.
#' Default is TRUE. Alternatively (movie=FALSE), you could use this function to get the node labels returned for each permutation starting with a perturbed variable.
#' @export mle.getPermMovie
#' @keywords probability
#' @keywords diffusion event
#' @keywords adaptive walk
#' @examples
#' mle.getPermMovie(patient.sig.nodes, patient, ig)
mle.getPermMovie = function(patient.sig.nodes, patient, ig, output_filepath, movie=TRUE) {
  # TODO :: Test on graphs. Do we like this design decision to subset neighborhood to get better view of network walk????
  ig.subset = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% names(unlist(neighborhood(ig, nodes=patient.sig.nodes, order=1))))])
  data.pvalsTmp = data.pvals[which(rownames(data.pvals) %in% V(ig.subset)$name),]
  V(ig.subset)$label.cex = 2
  V(ig.subset)$label = rep("", length(V(ig.subset)$name))
  tmp = get.adjacency(ig.subset, attr="weight")
  tmp = abs(tmp)
  igraphTestG2 = graph.adjacency(tmp, mode="undirected", weighted=TRUE)
  coords = layout.fruchterman.reingold(ig.subset)
  igraphObjectG = vector(mode="list", length=length(V(ig.subset)$name))
  names(igraphObjectG) = V(ig.subset)$name
  degs = list(degree(ig.subset))
  adjacency_matrix = list(as.matrix(get.adjacency(ig.subset, attr="weight")))

  # Phase 1: Do adaptive permutations ahead of time for each possible startNode in subGraphS.
  if (movie==TRUE) {
    V(ig.subset)$color = rep("white", length(igraphObjectG))
    V(ig.subset)$color[which(V(ig.subset)$name %in% patient.sig.nodes)] = "green"
    V(ig.subset)$label = rep("", length(igraphObjectG))
    current_node_set = NULL
    png(sprintf("%s/diffusionP1Movie%d_pt%d.png", output_filepath, length(current_node_set), patient), 500, 500)
    plot.igraph(ig.subset, layout=coords, vertex.color=V(ig.subset)$color,
                vertex.label=V(ig.subset)$label, edge.width=20*abs(E(ig.subset)$weight),
                vertex.size=round(10*(1-data.pvalsTmp[V(ig.subset)$name,patient]), 0), mark.col="dark red",
                mark.groups = names(current_node_set))
    title(sprintf("{%s}", paste(as.numeric(names(current_node_set) %in% patient.sig.nodes), collapse="")), cex.main=2)
    legend("bottomleft", legend=c("Jackpot Nodes", "Drawn Nodes"), fill=c("green", "dark red"))
    dev.off()
  }
  # Get all node names, so we know what all possible startNodes are.
  permutationByStartNode = list()
  bitStrings.pt = list()
  for (n in 1:length(patient.sig.nodes)) {
    current_node_set = mle.getPermN(n)
    print(sprintf("Calculating permutation %d of %d for patient %d...", n, length(patient.sig.nodes), patient))
    # Draw all nodes in graph
    current_node_set = NULL
    stopIterating=FALSE;
    startNode = patient.sig.nodes[n]
    hits = startNode
    numMisses = 0
    currentGraph = igraphObjectG

    while (stopIterating==FALSE) {
      current_node_set = c(current_node_set, startNode)
      #STEP 4: Diffuse p0 and p1 to connected nodes from current draw, startNode node
      #For unseen nodes, clear probabilities and add probability (p0/#unseen nodes)
      for (t in 1:length(currentGraph)) {
        currentGraph[[t]] = 0 #set probabilities of all nodes to 0
      }
      #determine base p0 probability
      baseP = p0/(length(currentGraph)-length(current_node_set))
      for (t in 1:length(currentGraph)) {
        if (!(names(currentGraph[t]) %in% current_node_set)) {
          currentGraph[[t]] = baseP  #set probabilities of unseen nodes to diffused p0 value, baseP
        } else {
          currentGraph[[t]] = 0
        }
      }
      p0_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
      currentGraph = graph.diffuseP1(p1, startNode, currentGraph, currentGraph[current_node_set], 1, verbose=FALSE)
      p1_event = sum(unlist(currentGraph[!(names(currentGraph) %in% current_node_set)]))
      if (abs(p1_event-1)>thresholdDiff) {
        extra.prob.to.diffuse = 1-p1_event
        currentGraph[names(current_node_set)] = 0
        currentGraph[!(names(currentGraph) %in% names(current_node_set))] = unlist(currentGraph[!(names(currentGraph) %in% names(current_node_set))]) + extra.prob.to.diffuse/sum(!(names(currentGraph) %in% names(current_node_set)))
      }
      #Set startNode to a node that is the max probability in the new currentGraph
      maxProb = names(which.max(currentGraph))

      if (movie==TRUE) {
        png(sprintf("%s/diffusionP1Movie%d_pt%d_%d-1.png", output_filepath, n, patient, length(current_node_set)), 500, 500)
        plot.igraph(ig.subset, layout=coords, vertex.color=V(ig.subset)$color,
                    vertex.label=V(ig.subset)$label, edge.width=20*abs(E(ig.subset)$weight),
                    vertex.size=5+round(50*sumHits, 0), mark.col="dark red",
                    mark.groups = names(current_node_set))
        title(sprintf("{%s}", paste(as.numeric(names(current_node_set) %in% patient.sig.nodes), collapse="")), cex.main=2)
        legend("bottomleft", legend=c("Jackpot Nodes", "Drawn Nodes"), fill=c("green", "dark red"))
        dev.off()
      }

      # Break ties: When there are ties, choose the first of the winners.
      startNode = names(currentGraph[maxProb[1]])
      if (length(c(startNode,current_node_set))>=(length(V(igraphTestG)$name))) {
        current_node_set = c(current_node_set, startNode)
        stopIterating = TRUE
      }

      if (movie==TRUE) {
        png(sprintf("%d/diffusionP1Movie%d_pt%d_%d-2.png", output_filepath, n, patient, length(current_node_set)), 500, 500)
        plot.igraph(ig.subset, layout=coords, vertex.color=V(ig.subset)$color,
                    vertex.label=V(ig.subset)$label, edge.width=20*abs(E(ig.subset)$weight),
                    vertex.size=5+round(50*sumHits, 0), mark.col="dark red",
                    mark.groups = c(names(startNode), names(current_node_set)))
        title(sprintf("{%s}", paste(as.numeric(c(names(current_node_set), names(startNode)) %in% patient.sig.nodes), collapse="")), cex.main=2)
        legend("bottomleft", legend=c("Jackpot Nodes", "Drawn Nodes"), fill=c("green", "dark red"))
        dev.off()
      }
    }

    permutationByStartNode[[n]] = current_node_set
    bitStrings.pt[[n]] = as.numeric(names(current_node_set) %in% patient.sig.nodes)
  }
  names(permutationByStartNode) = patient.sig.nodes
  optimalBitString[[patient]] = paste(bitStrings.pt[[which.max(unlist(lapply(lapply(bitStrings.pt, function(i) which(i==1)), function(i) sum(i))))]], collapse="")
  names(optimalBitString[[patient]]) = paste(patient.sig.nodes, collapse="/")

  if (movie==TRUE) {
    optBS = as.numeric(unlist(strsplit(optimalBitString[[patient]], split="")))
    opT = sum(optBS)
    n=which(lapply(bitStrings.pt, function(i) paste(i, collapse=""))==optimalBitString[[patient]])
    startNode = patient.sig.nodes[which(lapply(bitStrings.pt, function(i) paste(i, collapse=""))==optimalBitString[[patient]])]
    current_node_set = permutationByStartNode[[which(patient.sig.nodes==names(startNode)[1])]]
    png(sprintf("%s/diffusionMovie%d_pt%d_summary.png", output_filepath, n, patient), 500, 500)
    plot.igraph(ig.subset, layout=coords, vertex.color=V(ig.subset)$color,
                vertex.label=V(ig.subset)$label, edge.width=20*abs(E(ig.subset)$weight),
                vertex.size=5+round(50*sumHits, 0), mark.col="dark red",
                mark.groups = c(names(startNode), names(current_node_set)))
    title(sprintf("compressed{%s}\ndirectly encoded{log2(choose(%d-%d, %d-%d))}",
                  paste(as.numeric(names(current_node_set)[1:which(optBS==1)[opT]] %in% patient.sig.nodes), collapse=""),
                  length(igraphObjectG), opT, length(patient.sig.nodes), opT), cex.main=2)
    legend("bottomleft", legend=c("Jackpot Nodes", "Drawn Nodes"), fill=c("green", "dark red"))
    dev.off()
  }

  return(permutationByStartNode)
}


