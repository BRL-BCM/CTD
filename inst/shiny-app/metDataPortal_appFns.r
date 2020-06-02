load("/Users/lillian.rosa/OneDrive/MacFiles/9thCommitteeMeeting/Clinical_paper/data/sysdata.rda")
load(system.file(sprintf("shiny-app/disMod.RData"), package = "CTD"))
modelChoices <<- tolower(unique(sapply(list.files(system.file("ranks/ind_ranks",package = "CTD")),function(x) sub("[0-9]+-ranks.RData","",x))))
getData = function(input) {
  print("called getData()...")
  if (input$raworZscore == "Raw") {
    data = .GlobalEnv$data_raw
  } else if (input$raworZscore == "Normalized") {
    data = .GlobalEnv$data_norm
  } else if (input$raworZscore == "Zscored") {
    data = .GlobalEnv$data_zscore
  }
  pts = as.character(unlist(sapply(input$showThese, function(i) cohorts_coded[[i]])))
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
justify <- function(x, hjust="center", vjust="center", draw=TRUE){
  w <- sum(x$widths)
  h <- sum(x$heights)
  xj <- switch(hjust,
               center = 0.5,
               left = 0.5*w,
               right=unit(1,"npc") - 0.5*w)
  yj <- switch(vjust,
               center = 0.5,
               bottom = 0.5*h,
               top=unit(1,"npc") - 0.5*h)
  x$vp <- viewport(x=xj, y=yj)
  if(draw) grid.draw(x)
  return(x)
}

getPathwayMap = function(input) {
  #' Generate pathway map with patient data superimposed.
  #' @param Pathway.Name - The name of the pathway map you want to plot patient data on.
  #' @param PatientID - An identifier string associated with the patient.
  #' @param patient.zscore - A named vector of metabolites with corresponding z-scores.
  #' @param scalingFactor - Integer associated with increase in node size.
  #' @param outputFilePath - The directory in which you want to store image files.
  zscore.data = .GlobalEnv$data_zscore[,-c(1:8)]

  if (length(input$ptIDs)==0) {
    return(list(pmap = NULL, colorbar = NULL))
  } else {
    PatientID = input$ptIDs
    scalingFactor <<- input$scalingFactor
    print(scalingFactor)
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
    V(template.ig)$title = V(template.ig)$label
    if(length(nms) != 0){
      patient.zscore = patient.zscore[which(names(patient.zscore) %in% nms)]
      granularity = 2
      blues = colorRampPalette(c("blue", "white"))(granularity*ceiling(abs(min(na.omit(patient.zscore))))+1)
      reds = colorRampPalette(c("white", "red"))(granularity*ceiling(max(abs(na.omit(patient.zscore)))))
      redblue = c(blues, reds[2:length(reds)])

      for (i in 1:length(node.labels)) {
        if (node.labels[i] %in% nms) {
          if (!is.na(patient.zscore[node.labels[i]])) {
            V(template.ig)$size[i] = 10*scalingFactor*ceiling(abs(patient.zscore[node.labels[i]]))
            V(template.ig)$title[i] = sprintf("%s\nz-score = %.2f",node.labels[i],patient.zscore[node.labels[i]])
            V(template.ig)$color[i] = redblue[1+granularity*(ceiling(patient.zscore[node.labels[i]])-ceiling(min(na.omit(patient.zscore))))]
          } else {
            V(template.ig)$size[i] = 10
            V(template.ig)$color[i] = "#D3D3D3"
          }
        } else {
          V(template.ig)$size[i] = 10
          V(template.ig)$color[i] = "#D3D3D3"
        }
      }
    }

    # V(template.ig)$label = capitalize(tolower(V(template.ig)$label))
    wrap_strings = function(vector_of_strings,width){
      as.character(sapply(vector_of_strings, FUN=function(x){
        paste(strwrap(x, width=width), collapse="\n")
      }))
    }
    V(template.ig)$label = wrap_strings(V(template.ig)$label, 15)


    #visualize via visNetwork
    vis.ig=template.ig
    names(vertex_attr(vis.ig))[1] = "igraph_id"
    vertex_attr(vis.ig,"id") = V(vis.ig)$name
    vertex_attr(vis.ig)[['shape']] = replace(vertex_attr(vis.ig)[['shape']], which(node.types %in% c("Enzyme","Unknown")), "square")
    vertex_attr(vis.ig)[['shape']] = replace(vertex_attr(vis.ig)[['shape']], which(node.types %in% c("Metabolite","Minor+Metabolite","Intermediate")), "dot")
    vertex_attr(vis.ig)[['shape']] = replace(vertex_attr(vis.ig)[['shape']], which(node.types %in% c("Cofactor")), "triangle")
    vertex_attr(vis.ig)[['shape']] = replace(vertex_attr(vis.ig)[['shape']], which(node.types %in% c("Label")), "box")

    scaler= diff(range(V(vis.ig)$y))/800
    V(vis.ig)$x=(V(vis.ig)$x-range(V(vis.ig)$x)[1])/scaler
    V(vis.ig)$y=(V(vis.ig)$y-range(V(vis.ig)$y)[1])/scaler

    V(vis.ig)$label.family = "sans"
    #V(vis.ig)$label.cex = diff(range(V(vis.ig)$y))*700/diff(range(V(vis.ig)$y))*1.34/800
    V(vis.ig)$label.cex = 1.5
    E(vis.ig)$width = 5
    E(vis.ig)$color = "lightgrey"
    #V(vis.ig)$value = vertex_attr(vis.ig)[['size']]

    if(input$pathwayMapId == "All"){
      V(vis.ig)$label[node.types!="Label"]=" "
      vis.ig = delete.vertices(vis.ig, v=V(vis.ig)$name[which(node.types %in% c("Class","FinalPathway"))])
    }else{
      vis.ig = delete.vertices(vis.ig, v=V(vis.ig)$name[which(node.types %in% c("Class","Label","FinalPathway"))])
    }


    lnodes <- data.frame(label = c("Enzyme\nor Unknown",
                                   "Metabolite or\nMinor Metabolite\n or Intermediate",
                                   "Cofactor"),
                         font.size = 15,
                         shape = c( "square","dot","triangle"),
                         color = list(background="lightgrey",border = "grey"))

    df=data.frame(Node.Types=c("Enzyme or\nUnknown","Metabolite or\nMinor Metabolite or\nIntermediate","Cofactor"),
                  x=c(1,2,3),
                  y=c(1,2,3))

    pmap=visIgraph(vis.ig,idToLabel=F,physics = F) %>%
      visOptions( autoResize = F,
                  height = "800px",
                  highlightNearest = list(enabled = T, degree = 2, hover = F)) %>%
      visNodes(size=V(vis.ig)$size) %>%
      # visLegend(addNodes = lnodes, useGroups = FALSE, position = "left",width = 0.1, zoom = TRUE) %>%
      visInteraction(tooltipDelay = 0, hover=T)

    lvis = ggplot(df, aes(x=x, y=y, shape=Node.Types)) +
      geom_point(fill="lightgrey", color="darkgrey", size=12) +
      scale_shape_manual(values=c(24, 22, 21)) +
      theme(legend.position="bottom",
            legend.background = element_blank(),
            legend.key=element_blank(),

            #legend.key.size = unit(4, "shape"),
            legend.text=element_text(size=13),
            legend.title=element_blank())
    lvis <- cowplot::get_legend(lvis)
    lvis=justify(lvis,hjust = "right", vjust = "center",draw = FALSE)

  }

  # Get colorbar
  if(length(nms) != 0){
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
    leg=cowplot::get_legend(cb)
  }else{
    leg = textGrob("No metabolites (|Zscore|>2) found.\nShowing template pathway map.",gp=gpar(col="black", fontsize=14))
  }
  leg=justify(leg,hjust = "center", vjust = "center",draw = FALSE)

  # return(list(pmap = list(src=svg_filename, contentType = 'image/svg+xml'), colorbar = leg))
  return(list(pmap=pmap,colorbar = leg,shapeleg = lvis))
}

getPatientReport = function(input) {
  # Must display RAW, Anchor and Z-score values for all patients in input$ptIDs.
  # If in Miller2015 data, there are no raw and anchor values.
  zscore.data = .GlobalEnv$data_zscore
  tmp = rownames(zscore.data)
  zscore.data = apply(as.matrix(zscore.data[,which(colnames(zscore.data) %in% input$ptIDs)]), 1, function(i) mean(na.omit(i)))
  names(zscore.data) = tmp

  #norm.data = .GlobalEnv$data_norm
  #tmp = rownames(norm.data)
  #norm.data = apply(as.matrix(norm.data[,which(colnames(norm.data) %in% input$ptIDs)]), 1, function(i) mean(na.omit(i)))
  #names(norm.data) = tmp

  #raw.data = .GlobalEnv$data_raw
  #tmp = rownames(raw.data)
  #raw.data = apply(as.matrix(raw.data[,which(colnames(raw.data) %in% input$ptIDs)]), 1, function(i) mean(na.omit(i)))
  #names(raw.data) = tmp

  # MetaboliteName  RawIonIntensity Anchor(CMTRX.5 median value)  Zscore
  #data = data.frame(Metabolite=character(), Raw=numeric(), Anchor=numeric(), Zscore=numeric(), stringsAsFactors = FALSE)
  data = data.frame(Metabolite=character(), Zscore=numeric(), stringsAsFactors = FALSE)
  for (row in 1:length(zscore.data)) {
    data[row, "Metabolite"] = names(zscore.data)[row]
   # if (names(zscore.data)[row] %in% names(norm.data)) {
   #  data[row, "Anchor"] = round(norm.data[which(names(norm.data)==names(zscore.data)[row])], 2)
   #} else {
   #   data[row, "Anchor"] = NA
   #}
   #if (names(zscore.data)[row] %in% names(raw.data)) {
   #   data[row, "Raw"] = round(raw.data[which(names(raw.data)==names(zscore.data)[row])], 2)
   #} else {
   #   data[row, "Raw"] = NA
   #}
    data[row, "Zscore"] = round(zscore.data[row], 2)
  }
  data = data[order(abs(data$Zscore), decreasing = TRUE),]

  # Remove mets that were NA in raw, norm and zscore
  #ind0 = intersect(intersect(which(is.na(data[,"Raw"])), which(is.na(data[,"Anchor"]))), which(is.na(data[,"Zscore"])))
  ind0 = which(is.na(data[,"Zscore"]))
  ind1 = grep("x - ", rownames(data))
  # Next, Remove mets that were NA in raw, but not in Anchor. These will be displayed in separate table.
  # Note, these values were imputed and therefore should not be included in patient report, but should
  # be noted that these metabolites were normally found.

  # Find metabolites that were not detected but are normally detected:
  #      Metabolites with higher fill rate (>80%) should be detected.
  if (any(grep("^HEP", input$ptIDs))) {
    refs = .GlobalEnv$data_zscore[-grep("x - ", rownames(data_zscore)), grep("HEP-REF", colnames(.GlobalEnv$data_zscore))]
    ref.fil = Miller2015$`Times identifed in all 200 samples`/200
  } else {
    refs = .GlobalEnv$data_zscore[-grep("x - ", rownames(data_zscore)), grep("EDTA-REF", colnames(.GlobalEnv$data_zscore))]
    ref.fil = apply(refs, 1, function(i) sum(is.na(i))/length(i))
  }
  tmp = ref.fil[ind0]
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
  if (length(ind0)>0) { data = data[-c(ind0, ind1),] }
  print(dim(data))

  # Order by Fill Rate, then by abs(Zscore)
  missingMets = missingMets[order(missingMets[,"Reference Fill Rate"], decreasing = TRUE),]
  class(data[,"Zscore"]) = "numeric"
  data = data[order(abs(data[,"Zscore"]), decreasing = TRUE), ]
  names(data) = c("Metabolite", "Z-score")

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
    ggtitle(sprintf("%s (%s)", input$metSelect, input$metClass)) + labs(x="Z-score", y="Count")
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

getColumn <<- function(df,colname,rown="rowname"){
  vcol=df[,colname]
  if(rown == "rowname"){
    names(vcol) = rownames(df)
  }else{
    names(vcol) = df[,rown]
  }
  return(vcol)
}

rampred <<- function(numdata){
  #brks = quantile(numdata, probs = seq(.05, .95, .05), na.rm = TRUE)
  brks = quantile(seq(min(numdata),max(numdata),diff(range(numdata))/20), probs = seq(.00, 1, .05), na.rm = TRUE)
  clrs = round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
  clrs = rev(clrs)
  return(list(brks=brks,clrs=clrs))
}

getPrankDf=function(input){
  # get disease cohort p-value ranking dataframe
  pts = cohorts_coded[[input$diag_nw_Class]]
  load(system.file(sprintf("shiny-app/model/ptRanks_kmx%s.RData", input$kmx), package = "CTD"))
  #load(sprintf("/Users/lillian.rosa/Downloads/CTD/inst/shiny-app/model/ptRanks_kmx%s.RData", input$kmx))
  df.pranks = sapply(match(pts,names(pt_ranks)), function(x) getColumn(pt_ranks[[x]],"ctd","model"))
  colnames(df.pranks)=pts
  df.pranks = apply(df.pranks,c(1,2), function(x) as.numeric(sprintf("%.3e",x)))
  model.ind = match(input$diag_nw_Class,rownames(df.pranks))
  if (is.na(model.ind)) {model.ind = 1}
  return(list(df.pranks=df.pranks,
              model.ind=model.ind,
              np=length(pts)))
}

getPtPrankDf=function(input){
  load(system.file(sprintf("shiny-app/model/ptRanks_kmx%s.RData", input$kmx), package = "CTD"))
  ptID = input$pt_nw_ID
  pt.df.pranks = pt_ranks[[ptID]]
  rownames(pt.df.pranks)=pt.df.pranks[,"model"]
  pt.df.pranks = pt.df.pranks[,-which("model" %in% colnames(pt.df.pranks))]
  pt.df.pranks = apply(pt.df.pranks,c(1,2), function(x) as.numeric(sprintf("%.3e",x)))
  diag.ind = match(input$diag_nw_Class,rownames(pt.df.pranks))
  if(input$sigOnly) {
    brks=rampred(pt.df.pranks[pt.df.pranks<0.05])[["brks"]]
    clrs=rampred(pt.df.pranks[pt.df.pranks<0.05])[["clrs"]]
  }else{
    brks=rampred(pt.df.pranks)[["brks"]]
    clrs=rampred(pt.df.pranks)[["clrs"]]
  }

  # temp=datatable(pt.df.pranks,
  #                options = list(scrollX = TRUE,fixedColumns = TRUE, #dom = 't',
  #                               pageLength = 20,
  #                               lengthMenu = c(10, 15, 20)))
  #
  # for (n in 1:ncol(pt.df.pranks)) {
  #   temp=formatStyle(temp,colnames(pt.df.pranks)[n],
  #                    backgroundColor = styleInterval(rampred(pt.df.pranks[,n])[["brks"]],rampred(pt.df.pranks[,n])[["clrs"]]))
  # }
  # pt.df.pranks=temp
  pt.df.pranks=datatable(pt.df.pranks,
                         options = list(scrollX = TRUE,fixedColumns = TRUE, #dom = 't',
                                        pageLength = 20,
                                        lengthMenu = c(10, 15, 20))) %>%
    formatStyle(colnames(pt.df.pranks),
                backgroundColor = styleInterval(brks,clrs))
  return(pt.df.pranks)
}

getVlength <<- function(input){
  ptID=input$pt_nw_ID
  getDiag=names(unlist(sapply(cohorts_coded,function(x) which(ptID %in% x))))
  model=input$bgModel
  kmx=input$kmx
  if(model==getDiag){
    fold=which(cohorts_coded[[getDiag]]==ptID)
    # load latent-embedding, pruned network that is learnt from the rest of the patients diagnosed with the same disease.
    if (system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,fold), package='CTD') != ""){
      ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,fold), package='CTD'))[['ig_pruned']]
    }else{
      ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,1), package='CTD'))[['ig_pruned']]
    }
  }else{
    ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,1), package='CTD'))[['ig_pruned']]
  }
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name

  data_mx = as.matrix(.GlobalEnv$data_zscore[which(rownames(.GlobalEnv$data_zscore) %in% V(ig)$name), ])
  data_mx = suppressWarnings(apply(data_mx, c(1,2), as.numeric))

  zmets=data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID]
  zmets=names(zmets[abs(zmets)>2])

  if ( input$RangeChoice == "Top K perturbed metabolites only"){
    e = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% names(data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID][1:kmx]))])
    e = delete.vertices(e, V(e)[degree(e) == 0] )
  }else if(input$RangeChoice == "Abnormal Z-scored metabolites only, regardless of K"){
    e = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% zmets)])
    e = delete.vertices(e, V(e)[degree(e) == 0] )
  }else if(input$RangeChoice == "All Detected Metabolites"){
    e = ig
  }else{
    print("No Range Selected")
  }
  return(length(V(e)$name))
}

getPtResult=function(input){
  ptID=input$pt_nw_ID
  sprintf("network visualization selected patient: %s",ptID)
  getDiag=names(unlist(sapply(cohorts_coded,function(x) which(ptID %in% x))))
  model=input$bgModel
  kmx=input$kmx
  if(model==getDiag){
    fold=which(cohorts_coded[[getDiag]]==ptID)
    # load latent-embedding, pruned network that is learnt from the rest of the patients diagnosed with the same disease.
    if (system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,fold), package='CTD') != ""){
      ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,fold), package='CTD'))[['ig_pruned']]
    }else{
      ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,1), package='CTD'))[['ig_pruned']]
    }
  }else{
    ig = loadToEnv(system.file(sprintf('networks/ind_foldNets/bg_%s_ind_fold%d.RData',model,1), package='CTD'))[['ig_pruned']]
  }
  # get "ig" derived adjacency matrix
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name
  #adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))
  data_mx = as.matrix(.GlobalEnv$data_zscore[which(rownames(.GlobalEnv$data_zscore) %in% V(ig)$name), ])
  data_mx = suppressWarnings(apply(data_mx, c(1,2), as.numeric))

  # using single-node diffusion
  kmx = input$kmx
  S = data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID][1:kmx] # top kmx perturbed metabolites in ptID's profile
  ranks = loadToEnv(system.file(sprintf('ranks/ind_ranks/%s%d-ranks.RData',toupper(model), 1), package='CTD'))[["permutationByStartNode"]]
  ranks = lapply(ranks, tolower)
  ptBSbyK = singleNode.getPtBSbyK(names(S), ranks) # encode nodes
  res = mle.getEncodingLength(ptBSbyK, NULL, ptID, G) # get encoding length
  mets = unique(c(names(S), names(ptBSbyK[[which.max(res[,"d.score"])]]))) # best co-perturbed metabolite set is the most compressed subset of nodes
  p.mets=2^-(res[which.max(res[,"d.score"]),"d.score"]-log2(nrow(res))) # p value of this "modular perturbation"
  #print(mets)
  sprintf("p.mets = %.3e",p.mets)

  # generate igraph for disease-relevant metabolites of the selected patient
  zmets=data_mx[order(abs(data_mx[,ptID]), decreasing = TRUE),ptID]
  zmets.red=names(zmets[zmets>2])
  zmets.blue=names(zmets[zmets<(-2)])
  zmets=names(zmets[abs(zmets)>2])

  if ( input$RangeChoice == "Top K perturbed metabolites only"){
    e = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% mets)])
    e = delete.vertices(e, V(e)[degree(e) == 0] )
  }else if(input$RangeChoice == "Abnormal Z-scored metabolites only, regardless of K"){
    e = delete.vertices(ig, v=V(ig)$name[-which(V(ig)$name %in% zmets)])
    e = delete.vertices(e, V(e)[degree(e) == 0] )
  }else if(input$RangeChoice == "All Detected Metabolites"){
    e = ig
  }else{
    print("No Range Selected")
  }

  # assign groups and make ColourScale

  # group=c(sapply(zmets.red,function(x) x="Zscore(>2.0)"),
  #         sapply(zmets.blue,function(x) x="Zscore(<-2.0)"),
  #         sapply(setdiff(V(ig)$name,c(zmets.red,zmets.blue)), function(x) x="-2.0 < Zscore < 2.0"))
  #
  # ColourScale <- 'd3.scaleOrdinal().domain(["Zscore(>2.0)", "Zscore(<-2.0)", "-2.0 < Zscore < 2.0"]).range(["#990000", "#000066", "#d3d3d3"]);'

  node_ptMod = V(e)$name  %in% mets
  node_disMod = V(e)$name  %in% disMod[[model]]
  node_both = node_ptMod & node_disMod

  group=rep("Metabolites in Neither Modules",length(V(e)$name))
  for (l in 1:length(V(e)$name)) {
    if (node_both[l]) {
      group[l] = "Metabolites in Both Modules"
    } else if (node_ptMod[l]) {
      group[l] = "Metabolites only in Patient Module" #"55185D"
    } else if (node_disMod[l]) {
      group[l] = "Metabolites only in Disease Module"
    } else {
      group[l] = "Metabolites in Neither Modules"
    }
  }
  names(group)=V(e)$name

  ColourScale <- 'd3.scaleOrdinal().domain(["Metabolites only in Patient Module", "Metabolites only in Disease Module", "Metabolites in Both Modules","Metabolites in Neither Modules"]).range(["7554A3", "96C93C", "ECB602","#d3d3d3"]);'
  borderColor = rep("#d3d3d3",length(V(e)$name))
  borderColor[V(e)$name %in% zmets.blue] = "lightblue"
  borderColor[V(e)$name %in% zmets.red] = "red"
  # reds = intersect(V(e)$name[which(V(e)$name %in% names(S))], names(S[which(S>0)]))
  # blues = intersect(V(e)$name[which(V(e)$name %in% names(S))], names(S[which(S<0)]))
  #zmets.red=zmets.red[!zmets.red %in% reds]
  #zmets.blue=zmets.blue[!zmets.blue %in% blues]

  #generate networkd3

  net_p=igraph_to_networkD3(e)
  net_p$nodes$group=sapply(as.character(net_p$nodes$name),function(x) group[x])

  net_p$nodes$nodesize=sapply(as.character(net_p$nodes$name),function(x) abs(data_mx[x,ptID])^1.5)

  net_p$links$value=abs(net_p$links$value)*150
  linkColor_ptMod=net_p$nodes$name[net_p$links$source+1] %in% mets & net_p$nodes$name[net_p$links$target+1] %in% mets
  linkColor_disMod=net_p$nodes$name[net_p$links$source+1] %in% disMod[[model]] & net_p$nodes$name[net_p$links$target+1] %in% disMod[[model]]
  linkColor_both = linkColor_ptMod & linkColor_disMod
  linkColor = rep("lightgrey", length(linkColor_ptMod))
  for (l in 1:length(linkColor)) {
    if (linkColor_both[l]) {
      linkColor[l] = "ECB602"
    } else if (linkColor_ptMod[l]) {
      linkColor[l] = "7554A3"#"55185D"
    } else if (linkColor_disMod[l]) {
      linkColor[l] = "96C93C"
    } else {
      linkColor[l] = "lightgrey"
    }
  }
  net_p$links$color=linkColor

  ptNetwork=forceNetwork(Nodes = net_p$nodes, charge = -90, fontSize = 20, colourScale = JS(ColourScale),
                         Links = net_p$links,
                         linkColour = net_p$links$color,
                         #linkDistance = 100,
                         Nodesize = 'nodesize',
                         Source = 'source', Target = 'target',NodeID = 'name',Group = 'group',Value = "value",zoom = T,
                         opacity = 0.9,
                         legend = T)

  ptNetwork$x$nodes$border = borderColor

  ptNetwork <- htmlwidgets::onRender(ptNetwork, 'function(el, x) { d3.selectAll("circle").style("stroke", d => d.border); }')
  return(list(mets=mets,p.mets=p.mets,ptNetwork=ptNetwork))
}
