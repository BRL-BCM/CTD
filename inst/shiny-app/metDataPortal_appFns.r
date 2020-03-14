
comparePatientModPerts = function(input) {
  ig = graphs[[input$diagnosis]][["ig"]]
  igraphObjectG = vector(mode="list", length=length(V(ig)$name))
  names(igraphObjectG) = V(ig)$name
  .GlobalEnv$igraphObjectG = igraphObjectG
  adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))
  .GlobalEnv$adjacency_matrix = adjacency_matrix

  permutationByStartNode = graphs[[input$diagnosis]][["permutationByStartNode"]]
  diagnoses = graphs[[input$diagnosis]][["diagnoses"]]

  data = graphs[[input$diagnosis]][["data"]]
  data = data[which(rownames(data) %in% V(ig)$name),]
  data.pvals = apply(data, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
  data.pvals = t(data.pvals)

  pt.BSbyK = graphs[[input$diagnosis]][["pt.BSbyK"]]

  results = data.frame(patientID=character(), diagnosis=character(), d.score=numeric(), subsetSize=numeric(), bits=numeric(), bitstring=character(), stringsAsFactors = FALSE)
  r=1
  for (pt in 1:length(pt.BSbyK)) {
    res = mle.getEncodingLength(pt.BSbyK[[pt]], data.pvals, names(pt.BSbyK)[pt])
    results[r, ] = c(names(pt.BSbyK)[pt], diagnoses[pt], res[order(res[,"d.score"], decreasing = TRUE),][1,c("d.score", "subsetSize", "IS.alt", "optimalBS")])
    r=r+1
  }
  results[,"bits"] = round(results[,"bits"], 2)
  pvals = round(2^-results[,"d.score"], 4)
  pvals[pvals>1] = 1
  results[,"d.score"] = sprintf("%.2f (%.2f)", pvals, round(results[,"d.score"], 2))
  colnames(results) = c("patientID", "diagnosis", "p-value (d.score)", "# in Module", "# Bits", "Bitstring")

  return(results)
}

getMDS = function(input) {
  patientSim = graphs[[input$diagnosis]][["patientSim"]]
  d = graphs[[input$diagnosis]][["diagnoses"]]
  if (input$diagnosis=="zsd") {
    d = c(rep("PEX1", 19), rep("PEX7", 21), rep("negCntl", 40))
  }
  if (input$dim==2) {
    fitSim = cmdscale(patientSim, eig=FALSE, k=2)
    x = round(fitSim[,1], 2)
    y = round(fitSim[,2], 2)
    df = data.frame(x=x, y=y, color=d)
    mds_plot = plot_ly(df, x=~x, y=~y, color=~color, marker = list(size = 20))
  } else {
    fitSim = cmdscale(patientSim, eig=FALSE, k=3)
    x = round(fitSim[,1], 2)
    y = round(fitSim[,2], 2)
    z = round(fitSim[,3], 2)
    df = data.frame(x=x, y=y, z=z, color=d, label=colnames(patientSim))
    mds_plot = plot_ly(df, x=~x, y=~y, z=~z, color=~color, text=~label)
  }
  return(mds_plot)
}

extractModPerts = function(input) {
  .GlobalEnv$data = data(Miller2015)
  .GlobalEnv$ig = graphs[[input$diagnosis]][["ig_pruned"]]
  #.GlobalEnv$permutationByStartNode = graphs[[input$diagnosis]][["permutationByStartNode"]]
  #.GlobalEnv$diagnoses = graphs[[input$diagnosis]][["diagnoses"]]
  #.GlobalEnv$igraphObjectG = graphs[[input$diagnosis]][["igraphObjectG"]]
  #.GlobalEnv$adjacency_matrix = graphs[[input$diagnosis]][["adjacency_matrix"]]
  #.GlobalEnv$pt.BSbyK = graphs[[input$diagnosis]][["pt.BSbyK"]]

  data = data[which(rownames(data) %in% V(ig)$name),which(colnames(data) %in% c(input$ptID, input$ptID2))]
  data.pvals = apply(data, c(1,2), function(i) 2*pnorm(abs(i), lower.tail = FALSE))
  data.pvals = t(data.pvals)
  print(dim(data))

  simTbl = data.frame(kMets=numeric(), NCD=numeric(), DirSim=numeric(), Jaccard=numeric(), MutualInfoPer=numeric(), stringsAsFactors = FALSE)
  ptID = input$ptID
  kmax1 = length(pt.BSbyK[[ptID]])
  ptID2 = input$ptID2
  kmax2 = length(pt.BSbyK[[ptID2]])
  kmxx = min(kmax1, kmax2)-1
  print(kmxx)
  for (k in 1:kmxx) {
    tmp = mle.getPatientSimilarity(ig, pt.BSbyK[[ptID]][k], ptID, pt.BSbyK[[ptID2]][k], ptID2, data.pvals)
    simTbl[k,"kMets"] = k+1
    simTbl[k,"NCD"] = round(tmp$NCD, 2)
    simTbl[k,"DirSim"] = round(tmp$dirSim, 2)
    simTbl[k,"Jaccard"] = round(tmp$jacSim, 2)
    simTbl[k,"MutualInfoPer"] = round(tmp$mutualInfoPer, 2)
  }

  kmx = as.numeric(input$kmax)
  pt1.sig.nodes = names(which(pt.BSbyK[[input$ptID]][[kmx-1]]==1))
  pt2.sig.nodes = names(which(pt.BSbyK[[input$ptID2]][[kmx-1]]==1))
  sig.nodes = unique(c(pt1.sig.nodes, pt2.sig.nodes))
  print(sig.nodes)

  p.ig = induced_subgraph(ig, vids = sig.nodes)
  lnks = get.edgelist(p.ig, names=FALSE)
  lnks = as.data.frame(lnks)
  lnks = lnks - 1
  lnks = cbind(lnks, 10*abs(E(p.ig)$weight))
  colnames(lnks) = c("Source", "Target", "Weight")
  nds = as.data.frame(V(p.ig)$name)
  bm = rep("", length(V(p.ig)$name))
  bm[which(V(p.ig)$name %in% pt1.sig.nodes)] = "Patient1"
  bm[which(V(p.ig)$name %in% pt2.sig.nodes)] = "Patient2"
  bm[intersect(which(V(p.ig)$name %in% pt1.sig.nodes),
               which(V(p.ig)$name %in% pt2.sig.nodes))] = "Both"
  nds = cbind(nds, bm, rep(""))
  colnames(nds) = c("Identifier", "Group")
  ColourScale1 = 'd3.scaleOrdinal() .domain(["Patient1", "Patient2", "Both", ""]) .range(["#FF0000", "#00FF00", "#989848", "#FFFFFF"]);'
  pt_ig=forceNetwork(Links = lnks, Nodes = nds, Source = "Source", Target = "Target", Value = "Weight",
                     NodeID = "Identifier", Group = "Group", opacity = 1.0, zoom=TRUE, legend=TRUE, colourScale = JS(ColourScale1))

  return(list(sim=simTbl, pt_ig=pt_ig))
}






getData = function(input) {
  print("called getData()...")
  if (input$raworZscore == "Raw") {
    data = .GlobalEnv$all_raw_data
  } else if (input$raworZscore == "Normalized") {
    data = .GlobalEnv$all_norm_data
  } else if (input$raworZscore == "Zscored") {
    data = .GlobalEnv$all_data
  }
  pts = as.character(unlist(sapply(input$showThese, function(i) cohorts[[i]])))
  ind = which(colnames(data) %in% pts)
  data = data[,ind]
  data = apply(data, c(1,2), function(i) round(i, 2))
  res = cbind(rownames(data), data)
  colnames(res) = c("Metabolite", colnames(data))
  res = as.matrix(res)
  return(res)
}

getMetList = function(input) {
  # First, get rid of metabolites that have below fil rate
  data = Miller2015[,grep("IEM_", colnames(Miller2015))]
  ref = data[,which(diagnoses$diagnosis=="No biochemical genetic diagnosis")]
  ref.fil = apply(ref, 1, function(i) 1-(sum(is.na(i))/length(i)))
  ref = ref[which(ref.fil>0.66),]
  metClass = .GlobalEnv$metClass[which(ref.fil>0.66)]

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
getRefPop = function(input, norm.data) {
  print("getRefPop() called.")
  data = Miller2015[,grep("IEM_", colnames(Miller2015))]
  ref = data[,which(diagnoses$diagnosis=="No biochemical genetic diagnosis")]
  ref.fil = apply(ref, 1, function(i) 1-(sum(is.na(i))/length(i)))
  ref = ref[which(ref.fil>0.66),]
  print(dim(ref))
  print(input$metSelect)
  print(input$metSelect %in% rownames(ref))

  # Histogram plot of metabolite, with outliers that were removed during z-score calculation highlighted in red
  outlierSamples = which(ref[input$metSelect,] %in% boxplot.stats(ref[input$metSelect,])$out)
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

  per = list(up=length(which(.GlobalEnv$all_data[input$metSelect,]>2))/nrow(.GlobalEnv$all_data),
             down=length(which(neg(.GlobalEnv$all_data[input$metSelect,]) > 2))/nrow(.GlobalEnv$all_data))
  df = data.frame(Sample=1:ncol(.GlobalEnv$all_data), Zscore=.GlobalEnv$all_data[input$metSelect,])
  rare = ggplot(data=df, aes(x=Sample, y=Zscore)) + geom_point(size=1) +
    ggtitle(sprintf("Percentage with zscore >2 = %.2f.\nPercentage with zscore <-2 = %.2f.", per$up, per$down)) +
    geom_hline(yintercept=2, color="red") + geom_hline(yintercept=-2, color="red")

  x = ref[input$metSelect,-outlierSamples]
  d = qqnorm(as.numeric(x), plot.it = FALSE)
  xx = d$y
  zz = d$x
  t = lm(xx~zz, data=as.data.frame(x=xx, z=zz))
  mn.est = as.numeric(t$coefficients[1])
  sd.est = as.numeric(t$coefficients[2])

  samples = colnames(ref)[which(ref[input$metSelect,] %in% boxplot.stats(ref[input$metSelect,])$out)]
  values = ref[input$metSelect, which(ref[input$metSelect,] %in% boxplot.stats(ref[input$metSelect,])$out)]
  outlierSamples = cbind(samples, round(as.numeric(values),2))
  colnames(outlierSamples) = c("Samples Outliers", "Sample Value")

  return (list(hst=hst, outliers=outlierSamples, qq=qq, ests=list(mean=mn.est, std=sd.est), rare=rare, per=per))
}

getPatientReport = function(input) {
  # Must display RAW, Anchor and Z-score values for all patients in input$ptIDs. 
  # If in Miller2015 data, there are no raw and anchor values.
  zscore.data = Miller2015[,grep("IEM_", colnames(Miller2015))]
  ref = data[,which(diagnoses$diagnosis=="No biochemical genetic diagnosis")]
  ref.fil = apply(ref, 1, function(i) 1-(sum(is.na(i))/length(i)))

  # MetaboliteName  RawIonIntensity Anchor(CMTRX.5 median value)  Zscore
  data = data.frame(Metabolite=character(), Raw=numeric(), Anchor=numeric(), Zscore=numeric(), stringsAsFactors = FALSE)
  for (row in 1:length(raw.data)) {
    data[row, "Metabolite"] = names(raw.data)[row]
    #data[row, "Raw"] = round(raw.data[row], 2)
    #if (length(which(names(norm.data)==names(raw.data)[row]))>0) {
    #  data[row, "Anchor"] = round(norm.data[which(names(norm.data)==names(raw.data)[row])], 2)
    #} else {
    #  data[row, "Anchor"] = NA
    #}
    if (length(which(names(zscore.data)==names(raw.data)[row]))>0) {
      data[row, "Zscore"] = round(zscore.data[which(names(zscore.data)==names(raw.data)[row])], 2)
    } else {
      data[row, "Zscore"] = NA
    }
  }

  # Remove mets that were NA in raw, norm and zscore AND
  # Next, Remove mets that were NA in raw, but not in Anchor. These will be displayed in separate table.
  # Note, these values were imputed and therefore should not be included in patient report, but should
  # be noted that these metabolites were normally found.
  # Last, remove mets that were NA in raw, but not in Zscore.
  #ind0 = intersect(intersect(which(is.na(data[,"Raw"])), which(is.na(data[,"Anchor"]))), which(is.na(data[,"Zscore"])))
  #ind1 = intersect(which(is.na(data[,"Raw"])), which(!is.na(data[,"Anchor"])))
  #ind2 = intersect(which(is.na(data[,"Raw"])), which(!is.na(data[,"Zscore"])))
  #ind_all = data[,"Metabolite"][unique(c(ind0, ind1))] #ind2
  #tmp = ref.fil[ind_all]
  #report_these = tmp[which(tmp>0.80)]
  # Report these metabolites
  #missingMets = data.frame(Metabolite=character(), Reference.FillRate=numeric(), stringsAsFactors = FALSE)
  #if (length(report_these)>0) {
  #  for (i in 1:length(report_these)) {
  #    met = names(report_these)[i]
  #    missingMets[i, "Metabolite"] = met
  #    missingMets[i, "Reference.FillRate"] = ref.fil[which(names(ref.fil)==met)]
  #  }
  #  colnames(missingMets) = c("Compound", "Reference Fill Rate")
  #} else {
  #  missingMets = NULL
  #}
  #if (length(ind_all)>0) { data = data[-unique(c(ind0, ind1)),] }
  #print(dim(data))

  # Order by Fill Rate, then by abs(Zscore)
  #missingMets = missingMets[order(missingMets[,"Reference Fill Rate"], decreasing = TRUE),]
  #class(data[,"Zscore"]) = "numeric"
  #data = data[order(abs(data[,"Zscore"]), decreasing = TRUE), ]
  #names(data) = c("Metabolite", "Raw Ion Intensity", "Anchor", "Z-score")

  return(list(patientReport=data, missingMets=NULL))
}

getPathwayMap = function(input) {
  #' Generate pathway map with patient data superimposed.
  #' @param Pathway.Name - The name of the pathway map you want to plot patient data on.
  #' @param PatientID - An identifier string associated with the patient.
  #' @param patient.zscore - A named vector of metabolites with corresponding z-scores.
  #' @param scalingFactor - Integer associated with increase in node size.
  #' @param outputFilePath - The directory in which you want to store image files.
  zscore.data = Miller2015[, grep("IEM_", colnames(Miller2015))]
  
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
    res = plot.getPathwayIgraph(input, Pathway.Name)
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
      if ((node.labels[i] %in% nms)) {
        #if (!is.na(patient.zscore[node.labels[i]])) {
          #V(template.ig)$size[i] = scalingFactor*ceiling(abs(patient.zscore[node.labels[i]]))
          #V(template.ig)$color[i] = redblue[1+granularity*(ceiling(patient.zscore[node.labels[i]])-ceiling(min(na.omit(patient.zscore))))]
        #} else {
        #  V(template.ig)$size[i] = 1
        #  V(template.ig)$color[i] = "#D3D3D3"
        #}
      } else {
        print(i)
        V(template.ig)$size[i] = 1
        V(template.ig)$color[i] = "#D3D3D3"
      }
    }
    V(template.ig)$size[which(node.types=="Class")]
    V(template.ig)$label = capitalize(tolower(V(template.ig)$label))
    wrap_strings = function(vector_of_strings,width){
      as.character(sapply(vector_of_strings, FUN=function(x){
        paste(strwrap(x, width=width), collapse="\n")
      }))
    }
    V(template.ig)$label = wrap_strings(V(template.ig)$label, 15)
    V(template.ig)$label.cex = 0.75
    template.ig = delete.vertices(template.ig, v=grep(unlist(strsplit(Pathway.Name, split="-"))[1], V(template.ig)$label))

    svg_filename = system.file("shiny-app/metDataPortal_appFns.r", package="CTD")
    svg_filename = gsub("/metDataPortal_appFns.r", "", svg_filename)
    svg_filename = sprintf("%s/pmap-%s_%s.svg", svg_filename, Pathway.Name, input$diagClass)
    svg(filename = svg_filename, width=10, height=5)
    par(mar=c(1,0.2,1,1))
    plot.igraph(template.ig, layout=cbind(V(template.ig)$x, V(template.ig)$y), edge.arrow.size = 0.01, edge.width = 1,
                vertex.frame.color=V(template.ig)$color, main = gsub("-", " ", Pathway.Name))
    legend('bottom',legend=1:max(ceiling(V(template.ig)$size/scalingFactor)),
           pt.cex=seq(1, ceiling(max(V(template.ig)$size)), scalingFactor),
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





