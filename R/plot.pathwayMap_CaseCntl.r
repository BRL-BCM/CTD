#' Generate pathway map with patient perturbation data superimposed.
#'
#' @param Pathway - The name of the pathway map you want to plot patient data on.
#' @param groupName - An identifier string associated with the group being compared to a "baseline" condition.
#' @param pvalues - A named vector of metabolites with corresponding p-values (e.g., from a t-test). 
#' @param out.path - The directory in which you want to store image files.
#' @param SVG - Save as SVG or PNG? If SVG is TRUE, then an SVG image is saved. If FALSE, a PNG is saved.
#' @import Hmisc
#' @export plot.pathwayMap_CaseCntl
#' @usage plot.pathwayMap_CaseCntl(Pathway, groupName, pvalues, out.path, SVG=TRUE)
#' @examples
#' pwys = pathway.ListMaps_metabolon()
#' data(Miller2015)
#' Miller2015 = Miller2015[,grep("IEM", colnames(Miller2015))]
#' type = diagnoses$diagnosis
#' sdd = apply(Miller2015[,which(type=="Argininemia")], 1, sd)
#' Miller2015 = Miller2015[-which(sdd==0),]
#' pvals = unlist(apply(Miller2015, 1, function(i) t.test(i[which(type=="Argininemia")], i[which(type=="No biochemical genetic diagnosis")])$p.value))
#' names(pvals) = rownames(Miller2015)
#' plot.pathwayMap_CaseCntl("allPathways", groupName="Argininemia", pvalues=pvals, out.path=getwd(), SVG=TRUE)
plot.pathwayMap_CaseCntl = function(Pathway, groupName, pvalues, out.path, SVG=TRUE) {
  if (Pathway != "All Pathways") Pathway= paste0(unlist(strsplit(Pathway," ")),collapse = "-")  else Pathway ="allPathways"
  
  if (length(which(is.na(pvalues)))>0) {
    pvalues = pvalues[-which(is.na(pvalues))]
  }
  names(pvalues) = tolower(trimws(names(pvalues)))
  names(pvalues) = gsub("\\\"", "", names(pvalues))
  names(pvalues) = gsub("\\*", "", names(pvalues))
  load(system.file("extdata/complexNodes.RData", package = "CTD"))
  load(system.file(sprintf("extdata/RData/%s.RData", Pathway), package = "CTD"))
  template.g = ig

  if (length(which(names(pvalues)=="3-ureidopropionate"))>0) {
    names(pvalues)[which(names(pvalues)=="3-ureidopropionate")] = "ureidopropionate"
  }
  
  nodeDisplayNames= read.table(system.file(sprintf("extdata/%s/Name-%s.txt", Pathway, Pathway), package = "CTD"), header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = suppressWarnings(as.numeric(tmp.nms))
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeDisplayNames = as.character(tmp)
  names(nodeDisplayNames) = tmp.nms
  nodeDisplayNames = gsub("\\+", " ", nodeDisplayNames)
  # Load id to node types mappings
  nodeType = read.table(system.file(sprintf("extdata/%s/Type-%s.txt", Pathway, Pathway), package = "CTD"), header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = suppressWarnings(as.numeric(tmp.nms))
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeType = as.character(tmp)
  names(nodeType) = tmp.nms
  nodeType = nodeType[which(names(nodeType) %in% names(nodeDisplayNames))]

  node.labels = vector("character", length = length(V(template.g)$name))
  node.types = vector("character", length = length(V(template.g)$name))
  for (n in 1:length(V(template.g)$name)) {
    node.labels[n] = tolower(as.character(nodeDisplayNames[V(template.g)$name[n]]))
    node.types[n] = as.character(nodeType[V(template.g)$name[n]])
  }
  node.labels = as.character(sapply(node.labels, URLdecode))

  metabolon_to_data = read.csv(system.file("extdata/metabolon_to_data.txt", package = "CTD"), sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
  metabolon_to_data = metabolon_to_data[,c(1,2,3,4)]
  metabolon_to_data = apply(metabolon_to_data, 2, tolower)
  # Relabel nodes that have different names in dataset
  for (i in 1:nrow(metabolon_to_data)) {
    if (metabolon_to_data[i, "Data_Label"]!="") {
      lbl = metabolon_to_data[,"PathwayMap_Label"][i]
      ind = which(node.labels==lbl)
      if (length(ind)>0) {
        node.labels[ind] = metabolon_to_data[i, "Data_Label"]
      }
    }
  }

  nms = node.labels[which(node.labels %in% names(pvalues))]
  mapped = 1
  for (i in 1:length(node.labels)) {
    if (node.labels[i] %in% nms) {
      mapped = mapped + 1
      if (pvalues[which(names(pvalues)==node.labels[i])]<0.15) {
        # Yellow
        V(template.g)$size[i] = 2
        V(template.g)$color[i] = "#FFFF00"
        if (pvalues[which(names(pvalues)==node.labels[i])]< 0.05) {
          # Red
          V(template.g)$size[i] = 4
          V(template.g)$color[i] = "#FF0000"
        }
      } else {
        # Black
        V(template.g)$size[i] = 1
        V(template.g)$color[i] = "#000000"
      }
    } else {
      # Black
      V(template.g)$size[i] = 1
      V(template.g)$color[i] = "#000000"
    }
  }

  names(complexNodes) = tolower(names(complexNodes))
  mapped=0
  complexNodes = complexNodes[which(names(complexNodes) %in% node.labels)]
  # Next do the complex nodes
  if (length(which(names(complexNodes) %in% node.labels))>0) {
    for (n in 1:length(complexNodes)) {
      metsInComplex = as.character(sapply(complexNodes[[n]], tolower))
      metsInComplex = gsub("\\*", "", metsInComplex)
      mapped.mets = metsInComplex[which(metsInComplex %in% names(pvalues))]
      if (length(mapped.mets)>0) {
        print(sprintf("%s: %f", names(complexNodes)[n], length(mapped.mets)/length(metsInComplex)))
        mapped = mapped + length(mapped.mets)
        nodeSize = min(na.omit(pvalues[which(names(pvalues) %in% mapped.mets)]))
        if (is.na(nodeSize)) {
          # Black
          V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 1
          V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 1
          V(template.g)$color[which(node.labels==names(complexNodes[n]))] = "#000000"
        } else {
          if (nodeSize < 0.15) {
            V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 2
            V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 2
            V(template.g)$color[which(node.labels==names(complexNodes[n]))] = "#FFFF00"
            if (nodeSize < 0.05) {
              # Red
              V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 4
              V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 4
              V(template.g)$color[which(node.labels==names(complexNodes[n]))] = "#FF0000"
            }
          }
        }
      } else {
        V(template.g)$size[which(node.labels==names(complexNodes[n]))] = 1
        V(template.g)$size2[which(node.labels==names(complexNodes[n]))] = 1
        V(template.g)$color[which(node.labels==names(complexNodes[n]))] = "#000000"
      }
    }
  }

  V(template.g)$size[which(node.types=="Class")]
  V(template.g)$label = capitalize(tolower(V(template.g)$label))
  wrap_strings = function(vector_of_strings,width){
    as.character(sapply(vector_of_strings, FUN=function(x){
      paste(strwrap(x, width=width), collapse="\n")
    }))
  }
  V(template.g)$label = wrap_strings(V(template.g)$label, 15)
  V(template.g)$label.cex = 0.75
  template.g = delete.vertices(template.g, v=grep(unlist(strsplit(Pathway, split="-"))[1], V(template.g)$label))
  V(template.g)$color[which(V(template.g)$shape=="rectangle")] = rep("#32CD32", length(which(V(template.g)$shape=="rectangle")))

  if (Pathway=="allPathways") {
    V(template.g)$label.cex = 0.25
    V(template.g)$label = rep("", length(V(template.g)$name))
  }
  if (SVG) {
    svg(sprintf("%s/%s-%s.svg", out.path, Pathway, groupName), width=15, height=15)
    plot.igraph(template.g, layout=cbind(V(template.g)$x, V(template.g)$y), edge.arrow.size = 0.01, edge.width = 1,
                vertex.frame.color=V(template.g)$color, main = gsub("-", " ", Pathway))
    legend('bottom',legend=1:max(ceiling(V(template.g)$size)), col='black',pch=21, pt.bg='white', cex=2, horiz=TRUE)
    dev.off()
  } else {
    png(sprintf("%s/%s-%s.png", out.path, Pathway, groupName), width=10, height=10, units="in", res=300)
    plot.igraph(template.g, layout=cbind(V(template.g)$x, V(template.g)$y), edge.arrow.size = 0.01, edge.width = 1,
                vertex.frame.color=V(template.g)$color, main = gsub("-", " ", Pathway))
    legend('bottom',legend=1:max(ceiling(V(template.g)$size)), col='black',pch=21, pt.bg='white', cex=2, horiz=TRUE)
    dev.off()
  }

}


