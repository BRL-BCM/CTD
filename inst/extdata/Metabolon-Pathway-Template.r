# Make Metabolon Pathway template igraph

rm(list=ls())

allPwys = pathway.ListMaps_metabolon()
allPwys = allPwys[-1]
allPwys = as.character(sapply(allPwys, function(i) gsub(" ", "-", i)))

for (Pathway.Filename in allPwys) {
  # Load id to node display names mappings
  nodeDisplayNames= read.table(sprintf("/Users/lillian.rosa/Downloads/CTD/inst/extdata/%s/Name-%s.txt", 
                                       Pathway.Filename, Pathway.Filename), 
                               header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeDisplayNames, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = as.numeric(tmp.nms)
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeDisplayNames = as.character(tmp)
  names(nodeDisplayNames) = tmp.nms
  nodeDisplayNames = gsub("\\+", " ", nodeDisplayNames)
  # Load id to node types mappings
  nodeType = read.table(sprintf("/Users/lillian.rosa/Downloads/CTD/inst/extdata/%s/Type-%s.txt", Pathway.Filename, Pathway.Filename), 
                        header=TRUE, sep="\n", check.names = FALSE)
  tmp = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[2])
  tmp.nms = apply(nodeType, 1, function(i) unlist(strsplit(i, split= " = "))[1])
  ind = as.numeric(tmp.nms)
  ind2 = as.logical(sapply(ind, function(i) is.na(i)))
  tmp = tmp[-which(ind2)]
  tmp.nms = tmp.nms[-which(ind2)]
  nodeType = as.character(tmp)
  names(nodeType) = tmp.nms
  nodeType = nodeType[which(names(nodeType) %in% names(nodeDisplayNames))]
  
  # Crawl gml file and get NODE_NAME, x-coord, y-coord, width, height, fill color, type, outline color and outline-width
  df.metabolon = data.frame(id=numeric(),name=character(), x=numeric(), y=numeric(), size=numeric(), size2=numeric(),
                            color=character(), shape=character(), outline_color=character(), outline_width=numeric(), stringsAsFactors = FALSE)
  con=file(sprintf("/Users/lillian.rosa/Downloads/CTD/inst/extdata/%s/%s.gml", Pathway.Filename, Pathway.Filename), open="r")
  line=readLines(con)
  close(con)
  graphic.lns = grep("graphics",line)  
  edge.lns = grep("edge", line)
  graphic.lns = graphic.lns[which(graphic.lns < edge.lns[1])]
  # When it errors, it has finished loading node coords/ids and gotten to node edges
  for (l in graphic.lns) {
    id = unlist(strsplit(line[l-1], split="\t"))
    id = as.numeric(id[-which(id %in% c("", "id"))])
    
    x = unlist(strsplit(line[l+1], split="\t"))
    x = as.numeric(x[-which(x %in% c("", "x"))])
    
    y = unlist(strsplit(line[l+2], split="\t"))
    y = as.numeric(y[-which(y %in% c("", "y"))])
    
    w = unlist(strsplit(line[l+3], split="\t"))
    w = as.numeric(w[-which(w %in% c("", "w"))])
    
    h = unlist(strsplit(line[l+4], split="\t"))
    h = as.numeric(h[-which(h %in% c("", "h"))])
    
    fil = unlist(strsplit(line[l+5], split="\t"))
    fil = fil[-which(fil %in% c("", "fill"))]
    fil = gsub("\"", "", fil)
    
    sh = unlist(strsplit(line[l+6], split="\t"))
    sh = sh[-which(sh %in% c("", "type"))]
    sh = gsub("\"", "", sh)
    
    out_clr = unlist(strsplit(line[l+7], split="\t"))
    out_clr = out_clr[-which(out_clr %in% c("", "outline"))]
    out_clr = gsub("\"", "", out_clr)
    
    out_w = unlist(strsplit(line[l+8], split="\t"))
    out_w = as.numeric(out_w[-which(out_w %in% c("", "outline_width"))])
    
    name = unlist(strsplit(line[l+10], split="\t"))
    name = name[-which(name %in% c("", "label"))]
    name = gsub("\"", "", name)
    
    df.metabolon[which(graphic.lns==l),"id"] = id
    df.metabolon[which(graphic.lns==l),"name"] = name
    df.metabolon[which(graphic.lns==l),"x"] = x
    df.metabolon[which(graphic.lns==l),"y"] = y
    df.metabolon[which(graphic.lns==l),"size"] = w
    df.metabolon[which(graphic.lns==l),"size2"] = h
    df.metabolon[which(graphic.lns==l),"color"] = fil
    df.metabolon[which(graphic.lns==l),"shape"] = sh
    df.metabolon[which(graphic.lns==l),"outline_color"] = out_clr
    df.metabolon[which(graphic.lns==l),"outline_width"] = out_w
  }
  
  require(igraph)
  ig = make_empty_graph(n=0)
  for (node in 1:nrow(df.metabolon)) {
    ig = add_vertices(ig, nv=1, attr = as.list(df.metabolon[node,]))
  }
  # Get edges
  edge.lns = grep("edge", line)
  for (e in edge.lns) {
    source = unlist(strsplit(line[e+3], split="\t"))
    source = as.numeric(source[-which(source %in% c("", "source"))])
    
    target = unlist(strsplit(line[e+2], split="\t"))
    target = as.numeric(target[-which(target %in% c("", "target"))])
    
    ig = add_edges(ig, edges=c(V(ig)$name[which(V(ig)$id==source)], 
                               V(ig)$name[which(V(ig)$id==target)]))
  }
  for (n in 1:length(V(ig)$name)) {
    V(ig)$label[n] = as.character(nodeDisplayNames[V(ig)$name[n]])
  }
  table(V(ig)$shape)
  V(ig)$shape[which(V(ig)$shape=="ellipse")] = "circle"
  V(ig)$shape[which(V(ig)$shape=="hexagon")] = "rectangle"
  V(ig)$shape[which(V(ig)$shape=="octagon")] = "vrectangle"
  V(ig)$x = V(ig)$x + abs(min(V(ig)$x))
  V(ig)$y = V(ig)$y + abs(min(V(ig)$y))
  V(ig)$y = max(V(ig)$y) - V(ig)$y
  
  V(ig)$size = rep(1, length(V(ig)$name))
  V(ig)$size2 = rep(1, length(V(ig)$name))
  V(ig)$color = rep("grey", length(V(ig)$name))
  V(ig)$color[which(V(ig)$shape=="rectangle")] = rep("light green", length(which(V(ig)$shape=="rectangle")))
  V(ig)$label.cex = 0.5
  
  V(ig)$label = as.character(unlist(sapply(V(ig)$label, function(i) URLdecode(i))))
  #plot.igraph(ig, layout=cbind(V(ig)$x, V(ig)$y), edge.arrow.size = 0.01, edge.width = 1, vertex.frame.color=V(ig)$color)
  save(ig, file=sprintf("/Users/lillian.rosa/Downloads/CTD/inst/extdata/RData/%s.RData", Pathway.Filename))
}