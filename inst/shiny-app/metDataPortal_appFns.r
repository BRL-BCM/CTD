load("/Users/lillian.rosa/OneDrive/MacFiles/9thCommitteeMeeting/Clinical_paper/data/sysdata.rda")

getData = function(input) {
  print("called getData()...")
  if (input$raworZscore == "Raw") {
    data = .GlobalEnv$data_raw
  } else if (input$raworZscore == "Normalized") {
    data = .GlobalEnv$data_norm
  } else if (input$raworZscore == "Zscored") {
    data = .GlobalEnv$data_zscore
  }
  pts = as.character(unlist(sapply(input$showThese, function(i) cohorts[[i]])))
  ind = which(colnames(data) %in% pts)
  data = data[,ind]
  data = apply(data, c(1,2), function(i) round(i, 2))
  res = cbind(rownames(data), data)
  colnames(res) = c("Metabolite", colnames(data))
  res = as.matrix(res)
  
  print(sprintf("getData() outputted res dim = %d x %d", dim(res)[1], dim(res)[2]))
  
  return(res)
}

#### TAB 1 FUNCTIONS:  ####
getPathwayMap = function(input) {
  #' Generate pathway map with patient data superimposed.
  #' @param Pathway.Name - The name of the pathway map you want to plot patient data on.
  #' @param PatientID - An identifier string associated with the patient.
  #' @param patient.zscore - A named vector of metabolites with corresponding z-scores.
  #' @param scalingFactor - Integer associated with increase in node size.
  #' @param outputFilePath - The directory in which you want to store image files.
  zscore.data = .GlobalEnv$data_zscore[,-c(1:8)]
  
  if (length(input$ptIDs)==0) {
    return(list(pmap = list(src="", contentType = 'image/svg+xml'), colorbar = NULL))
  } else {
    PatientID = input$ptIDs
    scalingFactor = input$scalingFactor
    tmp = rownames(zscore.data)
    patient.zscore = as.matrix(zscore.data[,which(colnames(zscore.data) %in% input$ptIDs)])
    patient.zscore = apply(patient.zscore, 1, function(i) mean(na.omit(i)))
    names(patient.zscore) = tmp
    
    if (input$pathwayMapId=="All") { Pathway.Name = "allPathways" } else { Pathway.Name = gsub(" ", "-", input$pathwayMapId) } 
    res = shiny.getPathwayIgraph(input, Pathway.Name)
    template.ig = res$template.ig
    node.labels = res$nodeDisplayNames
    node.types = res$nodeType
    
    # Super-impose patient-specific (or mean cohort-specific) profiles onto metabolite nodes.
    # Ignore "complex nodes" (the squares)
    nms = node.labels[which(node.labels %in% names(patient.zscore))]
    patient.zscore = patient.zscore[which(names(patient.zscore) %in% nms)]
    granularity = 2
    blues = colorRampPalette(c("blue", "white"))(granularity*ceiling(abs(min(na.omit(patient.zscore))))+1)
    reds = colorRampPalette(c("white", "red"))(granularity*ceiling(max(abs(na.omit(patient.zscore)))))
    redblue = c(blues, reds[2:length(reds)])
    for (i in 1:length(node.labels)) {
      if (node.labels[i] %in% nms) {
        if (!is.na(patient.zscore[node.labels[i]])) {
          V(template.ig)$size[i] = scalingFactor*ceiling(abs(patient.zscore[node.labels[i]]))
          V(template.ig)$color[i] = redblue[1+granularity*(ceiling(patient.zscore[node.labels[i]])-ceiling(min(na.omit(patient.zscore))))]
        } else {
          V(template.ig)$size[i] = 1
          V(template.ig)$color[i] = "#D3D3D3"
        }
      } else {
        V(template.ig)$size[i] = 1
        V(template.ig)$color[i] = "#D3D3D3"
      }
    }
    V(template.ig)$label = capitalize(tolower(V(template.ig)$label))
    wrap_strings = function(vector_of_strings,width){
      as.character(sapply(vector_of_strings, FUN=function(x){
        paste(strwrap(x, width=width), collapse="\n")
      }))
    }
    V(template.ig)$label = wrap_strings(V(template.ig)$label, 15)
    V(template.ig)$label.cex = 0.75
    template.ig = delete.vertices(template.ig, v=V(template.ig)$name[which(V(template.ig)$shape %in% c("Label", "Class", "FinalPathway"))])
    
    svg_filename = system.file("shiny-app/metDataPortal_appFns.r", package="CTD")
    svg_filename = gsub("/metDataPortal_appFns.r", "", svg_filename)
    svg_filename = sprintf("%s/pmap-%s_%s.svg", svg_filename, Pathway.Name, input$diagClass)
    svg(filename = svg_filename, width=14, height=10)
    par(mar=c(1,0.2,1,1))
    plot.igraph(template.ig, layout=cbind(V(template.ig)$x, V(template.ig)$y), edge.arrow.size = 0.1, edge.width = 1,
                vertex.frame.color=V(template.ig)$color, main = gsub("-", " ", Pathway.Name))
    legend('bottom',legend=1:max(ceiling(V(template.ig)$size/scalingFactor)),
           pt.cex=seq(1, ceiling(max(V(template.ig)$size/scalingFactor)), 1),
           col='black',pch=21, pt.bg='white', cex=1, horiz=TRUE)
    dev.off()
    
    # Get colorbar
    z = seq(floor(min(na.omit(patient.zscore))), ceiling(max(na.omit(patient.zscore))), 1/granularity)
    df = data.frame(Zscores = z[1:length(redblue)],
                    Colors = redblue)
    if (length(which(apply(df, 1, function(i) any(is.na(i)))))>0) {
      df = df[-which(apply(df, 1, function(i) any(is.na(i)))),]
    }
    cb = ggplot(df, aes(x=1:nrow(df), y=Zscores, colour=Zscores)) + geom_point() + #ggtitle(input$diagClass) +
      scale_colour_gradient2(guide = "colourbar", low = "blue", mid="white", high="red") +
      guides(colour = guide_colourbar(draw.llim = min(df$Zscores), draw.ulim = max(df$Zscores),
                                      direction="horizontal", title.position = "top", barwidth = 10, barheight = 2, reverse = FALSE))
    g_legend=function(a.gplot){
      tmp = ggplot_gtable(ggplot_build(a.gplot))
      leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend = tmp$grobs[[leg]]
      return(legend)}
    leg = g_legend(cb);
    return(list(pmap = list(src=svg_filename, contentType = 'image/svg+xml'), colorbar = leg))
  }
}

getPatientReport = function(input) {
  # Must display RAW, Anchor and Z-score values for all patients in input$ptIDs. 
  # If in Miller2015 data, there are no raw and anchor values.
  zscore.data = .GlobalEnv$data_zscore
  tmp = rownames(zscore.data)
  zscore.data = apply(as.matrix(zscore.data[,which(colnames(zscore.data) %in% input$ptIDs)]), 1, function(i) mean(na.omit(i)))
  names(zscore.data) = tmp
  
  norm.data = .GlobalEnv$data_norm
  tmp = rownames(norm.data)
  norm.data = apply(as.matrix(norm.data[,which(colnames(norm.data) %in% input$ptIDs)]), 1, function(i) mean(na.omit(i)))
  names(norm.data) = tmp
  
  raw.data = .GlobalEnv$data_raw
  tmp = rownames(raw.data)
  raw.data = apply(as.matrix(raw.data[,which(colnames(raw.data) %in% input$ptIDs)]), 1, function(i) mean(na.omit(i)))
  names(raw.data) = tmp
  
  # MetaboliteName  RawIonIntensity Anchor(CMTRX.5 median value)  Zscore
  data = data.frame(Metabolite=character(), Raw=numeric(), Anchor=numeric(), Zscore=numeric(), stringsAsFactors = FALSE)
  for (row in 1:length(zscore.data)) {
    data[row, "Metabolite"] = names(zscore.data)[row]
    if (names(zscore.data)[row] %in% names(norm.data)) {
      data[row, "Anchor"] = round(norm.data[which(names(norm.data)==names(zscore.data)[row])], 2)
    } else {
      data[row, "Anchor"] = NA
    }
    if (names(zscore.data)[row] %in% names(raw.data)) {
      data[row, "Raw"] = round(raw.data[which(names(raw.data)==names(zscore.data)[row])], 2)
    } else {
      data[row, "Raw"] = NA
    }
    data[row, "Zscore"] = round(zscore.data[row], 2)
  }
  data = data[order(abs(data$Zscore), decreasing = TRUE),]
  
  # Remove mets that were NA in raw, norm and zscore
  ind0 = intersect(intersect(which(is.na(data[,"Raw"])), which(is.na(data[,"Anchor"]))), which(is.na(data[,"Zscore"])))
  # Next, Remove mets that were NA in raw, but not in Anchor. These will be displayed in separate table.
  # Note, these values were imputed and therefore should not be included in patient report, but should
  # be noted that these metabolites were normally found.

  # Find metabolites that were not detected but are normally detected
  ind1 = intersect(which(is.na(data[,"Anchor"])), which(!is.na(data[,"Zscore"])))
  if (any(grep("IEM_", input$ptIDs))) {
    refs = .GlobalEnv$data_zscore[, grep("HEP-REF", colnames(.GlobalEnv$data_zscore))]
  } else {
    refs = .GlobalEnv$data_zscore[, grep("EDTA-REF", colnames(.GlobalEnv$data_zscore))]
  }
  ref.fil = apply(refs, 1, function(i) sum(is.na(i))/length(i))
  tmp = ref.fil[ind1]
  report_these = tmp[which(tmp>0.80)]
  # Report these metabolites
  missingMets = data.frame(Metabolite=character(), Reference.FillRate=numeric(), stringsAsFactors = FALSE)
  if (length(report_these)>0) {
    for (i in 1:length(report_these)) {
      met = names(report_these)[i]
      missingMets[i, "Metabolite"] = met
      missingMets[i, "Reference.FillRate"] = ref.fil[which(names(ref.fil)==met)]
    }
    colnames(missingMets) = c("Compound", "Reference Fill Rate")
  } else {
    missingMets = NULL
  }
  if (length(ind_all)>0) { data = data[-ind0,] }
  print(dim(data))
  
  # Order by Fill Rate, then by abs(Zscore)
  missingMets = missingMets[order(missingMets[,"Reference Fill Rate"], decreasing = TRUE),]
  class(data[,"Zscore"]) = "numeric"
  data = data[order(abs(data[,"Zscore"]), decreasing = TRUE), ]
  names(data) = c("Metabolite", "Raw Ion Intensity", "Anchor", "Z-score")
  
  return(list(patientReport=data, missingMets=missingMets))
}



#### TAB 2 (INSPECT REFERENCE POPULATION) FUNCTIONS ####
getMetList = function(input) {
  # First, get rid of metabolites that have below fil rate
  ref.fil = 1-apply(.GlobalEnv$data_zscore[,-c(1:8)], 1, function(i) sum(is.na(i))/length(i))
  ref = data_zscore[which(ref.fil>0.66),]
  metClass = data_zscore[which(ref.fil>0.66), "SUPER_PATHWAY"]
  metClass[which(metClass=="")] = "Unknown"

  if (input$metClass=="Lipid") {
    return(rownames(ref)[which(metClass=="Lipid")])
  } else if (input$metClass=="Unknown") {
    return(rownames(ref)[which(metClass=="Unknown")])
  } else if (input$metClass=="Nucleotide") {
    return(rownames(ref)[which(metClass=="Nucleotide")])
  } else if (input$metClass=="Amino Acid") {
    return(rownames(ref)[which(metClass=="Amino Acid")])
  } else if (input$metClass=="Cofactors and Vitamins") {
    return(rownames(ref)[which(metClass=="Cofactors and Vitamins")])
  } else if (input$metClass=="Xenobiotics") {
    return(rownames(ref)[which(metClass=="Xenobiotics")])
  } else if (input$metClass=="Carbohydrate") {
    return(rownames(ref)[which(metClass=="Carbohydrate")])
  } else if (input$metClass=="Energy") {
    return(rownames(ref)[which(metClass=="Energy")])
  } else if (input$metClass=="Peptide") {
    return(rownames(ref)[which(metClass=="Peptide")])
  }
}

neg = function(x) { return(-x) }

# Get reference population statistics & plots
getRefPop = function(input) {
  print("getRefPop() called.")
  ref = .GlobalEnv$data_zscore[,-c(1:8)]
  if (input$anticoagulant=="EDTA") {
    ref = ref[,grep("EDTA-REF", colnames(ref))]
  } else {
    ref = ref[,grep("HEP-REF", colnames(ref))]
  }
  ref.fil = 1-apply(ref, 1, function(i) sum(is.na(i))/length(i))
  ref = ref[which(ref.fil>0.66),]
  print(dim(ref))
  print(input$metSelect)
  print(input$metSelect %in% rownames(ref))

  # Histogram plot of metabolite, with outliers that were removed during z-score calculation highlighted in red
  outlierSamples = which(as.numeric(ref[input$metSelect,]) %in% boxplot.stats(as.numeric(ref[input$metSelect,]))$out)
  print(sprintf("Length outlier samples = %d", length(outlierSamples)))
  if (length(outlierSamples)>0) {
    df = data.frame(x=as.numeric(ref[input$metSelect, -outlierSamples]))
  } else {
    df = data.frame(x=as.numeric(ref[input$metSelect,]))
  }
  hst = ggplot(data=df, aes(x=x)) + geom_histogram() +
    ggtitle(sprintf("%s (%s)", input$metSelect, input$metClass)) + labs(x="log(NormScaledImputed)", y="Count")
  y = quantile(df$x[!is.na(df$x)], c(0.05, 0.95))
  x = qnorm(c(0.05, 0.95))
  slope = diff(y)/diff(x)
  int = y[1L] - slope * x[1L]
  qq = ggplot(data=df, aes(sample=x)) + stat_qq() + ggtitle(sprintf("Normal QQ-Plot for %s (%s)", input$metSelect, input$metClass)) +
    geom_abline(slope = slope, intercept = int)

  per = list(up=length(which(ref[input$metSelect,]>2))/ncol(ref),
             down=length(which(neg(ref[input$metSelect,]) > 2))/ncol(ref))
  df = data.frame(Sample=1:ncol(ref), Zscore=as.numeric(ref[input$metSelect,]))
  rare = ggplot(data=df, aes(x=Sample, y=Zscore)) + geom_point(size=1) +
    ggtitle(sprintf("Percentage with zscore >2 = %.2f.\nPercentage with zscore <-2 = %.2f.", per$up, per$down)) +
    geom_hline(yintercept=2, color="red") + geom_hline(yintercept=-2, color="red")

  if (length(outlierSamples)>0) {x = ref[input$metSelect,-outlierSamples]}
  d = qqnorm(as.numeric(x), plot.it = FALSE)
  xx = d$y
  zz = d$x
  t = lm(xx~zz, data=as.data.frame(x=xx, z=zz))
  mn.est = as.numeric(t$coefficients[1])
  sd.est = as.numeric(t$coefficients[2])

  samples = colnames(ref)[which(as.numeric(ref[input$metSelect,]) %in% boxplot.stats(as.numeric(ref[input$metSelect,]))$out)]
  values = ref[input$metSelect, which(as.numeric(ref[input$metSelect,]) %in% boxplot.stats(as.numeric(ref[input$metSelect,]))$out)]
  outlierSamples = cbind(samples, round(as.numeric(values),2))
  colnames(outlierSamples) = c("Samples Outliers", "Sample Value")

  return (list(hst=hst, outliers=outlierSamples, qq=qq, ests=list(mean=mn.est, std=sd.est), rare=rare, per=per))
}




#### TAB 3 (NETWORK-ASSISTED DIAGNOSTICS) FUNCTIONS ####
getPrankDf=function(input){
  ptID=input$pt_nw_ID
  data_mx = .GlobalEnv$data_zscore[,-c(1:8)]
  
  df.pranks = data.frame(Disease_Model=character(),P_Value=character(),stringsAsFactors = FALSE)
  r=1
  S.pt = data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID]
  for(model in names(cohorts_coded)){
    ig=loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%s.RData',model,1), package='CTD'))[['ig_pruned']]
    G = vector(mode="list", length=length(V(ig)$name))
    names(G) = V(ig)$name
    adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))
    data_mx = data_mx[which(rownames(data_mx) %in% V(ig)$name), ]
    # load node ranks
    ranks = loadToEnv(system.file(sprintf('ranks/ind_ranks/%s%d-ranks.RData',toupper(model), 1), package='CTD'))[["permutationByStartNode"]]
    ranks = lapply(ranks, tolower)
    diag = input$diag_nw_Class
    S = S.pt[1:input$kmx]
    S = S[which(names(S) %in% names(G))]
    ptBSbyK = singleNode.getPtBSbyK(names(S), ranks, num.misses = log2(length(G)))
    res = mle.getEncodingLength(ptBSbyK, NULL, ptID, G)
    df.pranks[r,"Disease_Model"] = model
    df.pranks[r,"P_Value"] = sprintf("%.3e",2^-(res[which.max(res[,"d.score"]),"d.score"]-log2(nrow(res))))
    r=r+1
  }
  df.pranks=df.pranks[order(as.numeric(df.pranks[,"P_Value"])),]
  return(df.pranks)
}

getPtResult=function(input){
  ptID=input$pt_nw_ID
  print(ptID)
  getDiag=sapply(cohorts_coded,function(x) which(x==ptID))
  model=input$bgModel
  if(model==names(getDiag[sapply(getDiag,length)>0])){
    fold=getDiag[sapply(getDiag,length)>0] 
    # load latent-embedding, pruned network that is learnt from the rest of the patients diagnosed with the same disease.
    ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%s.RData',model,fold), package='CTD'))[['ig_pruned']]
  }else{
    ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%s.RData',model,1), package='CTD'))[['ig_pruned']]
  }
  # get "ig" derived adjacency matrix
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name
  adjacency_matrix <<- list(as.matrix(get.adjacency(ig, attr="weight")))
  data_mx = data_mx.og[which(rownames(data_mx.og) %in% V(ig)$name), ]
  # using single-node diffusion
  kmx = 30
  S = data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID][1:kmx] # top kmx perturbed metabolites in ptID's profile
  ranks = loadToEnv(system.file(sprintf('ranks/ind_ranks/%s%d-ranks.RData',toupper(model), 1), package='CTD'))[["permutationByStartNode"]]
  ranks = lapply(ranks, tolower)
  ptBSbyK = singleNode.getPtBSbyK(names(S), ranks) # encode nodes
  res = mle.getEncodingLength(ptBSbyK, NULL, ptID, G) # get encoding length
  mets = unique(c(names(S), names(ptBSbyK[[which.max(res[,"d.score"])]]))) # best co-perturbed metabolite set is the most compressed subset of nodes
  p.mets=2^-(res[which.max(res[,"d.score"]),"d.score"]-log2(nrow(res))) # p value of this "modular perturbation"
  print(mets)
  print(p.mets)
  
  # generate igraph for disease-relevant metabolites of the selected patient
  zmets=data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID]
  zmets.red=names(zmets[zmets>2])
  zmets.blue=names(zmets[zmets<(-2)])
  zmets=names(zmets[abs(zmets)>2])
  #e = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% zmets)])
  #e = delete.vertices(e, V(e)[degree(e) == 0] )
  #e = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% mets)])
  reds = intersect(V(e)$name[which(V(e)$name %in% names(S))], names(S[which(S>0)]))
  blues = intersect(V(e)$name[which(V(e)$name %in% names(S))], names(S[which(S<0)]))
  zmets.red=zmets.red[!zmets.red %in% reds]
  zmets.blue=zmets.blue[!zmets.blue %in% blues]
  #generate networkd3
  group=c(sapply(reds,function(x) x="pos_sig"),
          sapply(blues,function(x) x="neg_sig"),
          sapply(zmets.red,function(x) x="pos_nonsig"),
          sapply(zmets.blue,function(x) x="neg_nonsig"))
  net_p=igraph_to_networkD3(e)
  net_p$nodes$group=sapply(as.character(net_p$nodes$name),function(x) group[x])
  ColourScale <- 'd3.scaleOrdinal().domain(["pos_sig", "neg_sig","pos_nonsig","neg_nonsig"]).range(["#990000", "#000066","pink","#CCCCFF"]);'
  net_p$nodes$nodesize=sapply(as.character(net_p$nodes$name),function(x) abs(data_mx[x,ptID])^2)
  net_p$links$value=abs(net_p$links$value)*100
  linkColor=net_p$nodes$name[net_p$links$source+1] %in% mets & net_p$nodes$name[net_p$links$target+1] %in% mets
  net_p$links$color=ifelse(linkColor,"red","lightgrey")
  
  ptNetwork=forceNetwork(Nodes = net_p$nodes, charge = -90, fontSize = 20, colourScale = JS(ColourScale), 
                         Links = net_p$links,
                         linkColour = net_p$links$color,
                         #linkDistance = 100,
                         Nodesize = 'nodesize',
                         Source = 'source', Target = 'target',NodeID = 'name',Group = 'group',Value = "value",zoom = T,
                         opacity = 0.9,
                         legend = T)
  
  return(list(mets=mets,p.mets=p.mets,ptNetwork=ptNetwork))
}



