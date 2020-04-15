#' Metabolite set enrichment analysis (MSEA) using pathway knowledge curated by Metabolon
#'
#' A function that returns the pathway enrichment score for all perturbed metabolites in a patient's full metabolomic profile.
#' # Main MSEA Analysis Function that implements the entire methodology
#' This is a methodology for the analysis of global molecular profiles called Metabolite Set Enrichment Analysis (MSEA). It determines
#' whether an a priori defined set of metabolites shows statistically significant, concordant differences between two biological
#' states (e.g. phenotypes). MSEA operates on all metabolites from an experiment, rank ordered by the signal to noise ratio and
#' determines whether members of an a priori defined metabolite set are nonrandomly distributed towards the top or bottom of the
#' list and thus may correspond to an important biological process. To assess significance the program uses an empirical
#' permutation procedure to test deviation from random that preserves correlations between metabolites. For details see 
#' Subramanian et al 2005.
#' 
#' @param ds: Input metabolite expression dataset
#' @param cls:  Input class vector (phenotype)
#' @param met.db: Metabolite set database in GMT format
#' @param nperm: Number of random permutations (default: 1000)
#' @param weighted.score.type: Enrichment correlation-based weighting: 0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1)
#' @param adjust.FDR.q.val: Adjust the FDR q-vals (default: F)
#' @param met.size.threshold.min: Minimum size (in metabolites) for database metabolite sets to be considered (default: 25)
#' @param met.size.threshold.max: Maximum size (in metabolites) for database metabolite sets to be considered (default: 500)
#' @param random.seed: Random number metaboliterator seed. (default: 123456)
#' @return report: Global output report, sorted by NES in decreasing order.
#' @export shiny.getMSEA_Metabolon
#' @examples
#' # If the .GMT file isn't already created, create it.
#' data(Miller2015)
#' population = rownames(Miller2015)
#' paths.hsa = list.dirs(path=system.file("extdata", package="CTD"), full.names = FALSE)
#' paths.hsa = paths.hsa[-which(paths.hsa %in% c("", "RData", "allPathways", "MSEA_Datasets"))]
#' sink(sprintf("%s/Metabolon.gmt", system.file("extdata/MSEA_Datasets", package="CTD")))
#' for (p in 1:length(paths.hsa)) {
#'   load(system.file(sprintf("/extdata/RData/%s.RData", paths.hsa[p]), package="CTD"))
#'   pathway.compounds = V(ig)$label[which(V(ig)$shape=="circle")]
#'   pathCompIDs = unique(tolower(pathway.compounds[which(pathway.compounds %in% population)]))
#'   print(sprintf("%s         %s", paths.hsa[p], paste(pathCompIDs, collapse="    ")), quote=FALSE)
#' }
#' sink()
#' res = shiny.getMSEA_Metabolon(input, cohorts)
#'     # The format (columns) for the global result files is as follows.
#'     Pathway : Pathway name.
#'     SIZE : Size of the set in metabolites.
#'     NES : Normalized (multiplicative rescaling) normalized enrichment score.
#'     NOM p-val : Nominal p-value (from the null distribution of the metabolite set).
#'     FDR q-val: False discovery rate q-values
#'     FWER p-val: Family wise error rate p-values.
#'     glob.p.val: P-value using a global statistic (number of sets above the set's NES).
shiny.getMSEA_Metabolon = function(input, cohorts) {
  # First, get the class labels
  class.labels = colnames(.GlobalEnv$data_zscore[,-c(1:8)])
  class.labels[which(!(class.labels %in% cohorts[[input$diagClass]]))] = 0
  class.labels[which(class.labels %in% cohorts[[input$diagClass]])] = 1
  class.labels = as.numeric(class.labels)
  phen1 = 0
  phen2 = 1
  # Second, get the metabolomics profiling data, with metabolites as rows, and samples as columns.
  data = .GlobalEnv$data_zscore[,-c(1:8)]
  fill.rate = apply(data, 1, function(i) sum(is.na(i))/length(i))
  data = data[which(fill.rate<0.33), ]
  for (r in 1:nrow(data)) {
    data[r, which(is.na(data[r, ]))] = min(na.omit(as.numeric(data[r,])))
  }
  any(is.na(data))
  any(is.infinite(unlist(data)))
  rownames(data) = tolower(rownames(data))
  data = apply(data, c(1, 2), as.numeric)
  dim(data)
  
  # Reorder both data columns and class labels so diseased are first, then controls
  data = data[, order(class.labels, decreasing = TRUE)]
  class.labels = class.labels[order(class.labels, decreasing = TRUE)]
  # Third, provide the file extension local to the installation of the CTD package for the 
  # desired pathway knowledgebase .GMT.
  met.db = system.file("extdata/MSEA_Datasets/Metabolon.gmt", package="CTD")

  MSEA.MetaboliteRanking = function(A, class.labels, metabolite.labels, nperm, sigma.correction = "MetaboliteCluster", replace=F) {
    A = A + 0.00000001
    N = length(A[,1])
    Ns = length(A[1,])
    
    subset.mask = matrix(0, nrow=Ns, ncol=nperm)
    reshuffled.class.labels1 = matrix(0, nrow=Ns, ncol=nperm)
    reshuffled.class.labels2 = matrix(0, nrow=Ns, ncol=nperm)
    class.labels1 = matrix(0, nrow=Ns, ncol=nperm)
    class.labels2 = matrix(0, nrow=Ns, ncol=nperm)
    
    order.matrix = matrix(0, nrow = N, ncol = nperm)
    obs.order.matrix = matrix(0, nrow = N, ncol = nperm)
    s2n.matrix = matrix(0, nrow = N, ncol = nperm)
    obs.s2n.matrix = matrix(0, nrow = N, ncol = nperm)
    
    obs.metabolite.labels = vector(length = N, mode="character")
    obs.metabolite.descs = vector(length = N, mode="character")
    obs.metabolite.symbols = vector(length = N, mode="character")
    
    M1 = matrix(0, nrow = N, ncol = nperm)
    M2 = matrix(0, nrow = N, ncol = nperm)
    S1 = matrix(0, nrow = N, ncol = nperm)
    S2 = matrix(0, nrow = N, ncol = nperm)
    
    C = split(class.labels, class.labels)
    class1.size = length(C[[1]])
    class2.size = length(C[[2]])
    class1.index = seq(1, class1.size, 1)
    class2.index = seq(class1.size + 1, class1.size + class2.size, 1)
    for (r in 1:nperm) {
      class1.subset = sample(class1.index, size = ceiling(class1.size), replace = replace)
      class2.subset = sample(class2.index, size = ceiling(class2.size), replace = replace)
      class1.subset.size = length(class1.subset)
      class2.subset.size = length(class2.subset)
      subset.class1 = rep(0, class1.size)
      for (i in 1:class1.size) {
        if (is.element(class1.index[i], class1.subset)) {
          subset.class1[i] = 1
        }
      }
      subset.class2 = rep(0, class2.size)
      for (i in 1:class2.size) {
        if (is.element(class2.index[i], class2.subset)) {
          subset.class2[i] = 1
        }
      }
      subset.mask[, r] = as.numeric(c(subset.class1, subset.class2))
      fraction.class1 = class1.size/Ns
      fraction.class2 = class2.size/Ns
      # Random (unbalanced permutation)
      full.subset = c(class1.subset, class2.subset)
      label1.subset = sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] = rep(0, Ns)
      reshuffled.class.labels2[, r] = rep(0, Ns)
      class.labels1[, r] = rep(0, Ns)
      class.labels2[, r] = rep(0, Ns)
      for (i in 1:Ns) {
        m1 = sum(!is.na(match(label1.subset, i)))
        m2 = sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] = m1
        reshuffled.class.labels2[i, r] = m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] = m2
          class.labels2[i, r] = 0
        } else {
          class.labels1[i, r] = 0
          class.labels2[i, r] = m2
        }
      }
    }
    # compute S2N for the random permutation matrix
    P = reshuffled.class.labels1 * subset.mask
    n1 = sum(P[,1])
    M1 = A %*% P
    M1 = M1/n1
    A2 = A*A
    S1 = A2 %*% P
    S1 = S1/n1 - M1*M1
    if (n1>1) {
      S1 = sqrt(abs((n1/(n1-1)) * S1))
    }
    P = reshuffled.class.labels2 * subset.mask
    n2 = sum(P[,1])
    M2 = A %*% P
    M2 = M2/n2
    A2 = A*A
    S2 = A2 %*% P
    S2 = S2/n2 - M2*M2
    if (n2>1) {
      S2 = sqrt(abs((n2/(n2-1)) * S2))
    }
    if (sigma.correction == "MetaboliteCluster") {  # small sigma "fix" as used in MetaboliteCluster
      S2 = ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 = ifelse(S2 == 0, 0.2, S2)
      S1 = ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 = ifelse(S1 == 0, 0.2, S1)
    }
    M1 = M1 - M2
    S1 = S1 + S2
    s2n.matrix = M1/S1
    for (r in 1:nperm) {order.matrix[, r] = order(s2n.matrix[, r], decreasing=T)}
    # compute S2N for the "observed" permutation matrix
    P = class.labels1 * subset.mask
    n1 = sum(P[,1])
    if (n1>1) {
      M1 = A %*% P
      M1 = M1/n1
      A2 = A*A
      S1 = A2 %*% P
      S1 = S1/n1 - M1*M1
      S1 = sqrt(abs((n1/(n1-1)) * S1))
    } else {
      M1 = A %*% P
      A2 = A*A
      S1 = A2 %*% P
      S1 = S1 - M1*M1
      S1 = sqrt(abs(S1))
    }
    P = class.labels2 * subset.mask
    n2 = sum(P[,1])
    if (n2>1) {
      M2 = A %*% P
      M2 = M2/n2
      A2 = A*A
      S2 = A2 %*% P
      S2 = S2/n2 - M2*M2
      S2 = sqrt(abs((n2/(n2-1)) * S2))
    } else {
      M2 = A %*% P
      A2 = A*A
      S2 = A2 %*% P
      S2 = S2 - M2*M2
      S2 = sqrt(abs(S2))
    }
    if (sigma.correction == "MetaboliteCluster") {  # small sigma "fix" as used in MetaboliteCluster
      S2 = ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 = ifelse(S2 == 0, 0.2, S2)
      S1 = ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 = ifelse(S1 == 0, 0.2, S1)
    }
    M1 = M1 - M2
    S1 = S1 + S2
    obs.s2n.matrix = M1/S1
    for (r in 1:nperm) {obs.order.matrix[,r] = order(obs.s2n.matrix[,r], decreasing=T)            }
    
    return(list(s2n.matrix = s2n.matrix, obs.s2n.matrix = obs.s2n.matrix, order.matrix = order.matrix, obs.order.matrix = obs.order.matrix))
  }
  
  MSEA.EnrichmentScore2 = function(metabolite.list, metabolite.set, weighted.score.type = 1, correl.vector = NULL) {
    # Computes the weighted MSEA score of metabolite.set in metabolite.list. It is the same calculation as in
    # MSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
    # This call is intended to be used to asses the enrichment of random permutations rather than the
    # observed one.
    N = length(metabolite.list)
    Nh = length(metabolite.set)
    Nm =  N - Nh
    
    loc.vector = vector(length=N, mode="numeric")
    peak.res.vector = vector(length=Nh, mode="numeric")
    valley.res.vector = vector(length=Nh, mode="numeric")
    tag.correl.vector = vector(length=Nh, mode="numeric")
    tag.loc.vector = vector(length=Nh, mode="numeric")
    tag.diff.vector = vector(length=Nh, mode="numeric")
    
    loc.vector[metabolite.list] = seq(1, N)
    tag.loc.vector = loc.vector[metabolite.set]
    tag.loc.vector = sort(tag.loc.vector, decreasing = F)
    
    if (weighted.score.type == 0) {
      tag.correl.vector = rep(1, Nh)
    } else if (weighted.score.type == 1) {
      tag.correl.vector = correl.vector[tag.loc.vector]
      tag.correl.vector = abs(tag.correl.vector)
    } else if (weighted.score.type == 2) {
      tag.correl.vector = correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector = abs(tag.correl.vector)
    } else {
      tag.correl.vector = correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector = abs(tag.correl.vector)
    }
    tag.correl.vector[is.na(tag.correl.vector)] = 1
    
    norm.tag = 1.0/sum(tag.correl.vector)
    tag.correl.vector = tag.correl.vector * norm.tag
    norm.no.tag = 1.0/Nm
    
    tag.diff.vector[1] = (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] = tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector = tag.diff.vector * norm.no.tag
    peak.res.vector = cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector = peak.res.vector - tag.correl.vector
    max.ES = max(peak.res.vector)
    min.ES = min(valley.res.vector)
    ES = signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
    
    return(list(ES = ES))
  }
  
  print(" *** Running MSEA Analysis...")
  nperm=1000
  weighted.score.type=1
  adjust.FDR.q.val = F
  met.size.threshold.min = 5
  met.size.threshold.max = 500
  random.seed = 760435
  nom.p.val.threshold = -1
  fwer.p.val.threshold = -1
  fdr.q.val.threshold = -1
  set.seed(seed=random.seed, kind = NULL)
  adjust.param = 0.5
  metabolite.labels = rownames(data)
  sample.names = colnames(data)
  A = data.matrix(data)
  cols = ncol(A)
  rows = nrow(A)

  # Read input metabolite set database
  temp = readLines(met.db)
  max.Ng = length(temp)
  temp.size.G = vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) { temp.size.G[i] = length(unlist(strsplit(temp[[i]], "    "))) - 2 }
  
  max.size.G = max(temp.size.G)
  gs = matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
  temp.names = vector(length = max.Ng, mode = "character")
  temp.desc = vector(length = max.Ng, mode = "character")
  met.count = 1
  for (i in 1:max.Ng) {
    metabolite.set.size = length(unlist(strsplit(temp[[i]], "    "))) - 2
    met.line = noquote(unlist(strsplit(temp[[i]], "    ")))
    metabolite.set.name = met.line[1]
    metabolite.set.desc = met.line[2]
    metabolite.set.tags = vector(length = metabolite.set.size, mode = "character")
    for (j in 1:metabolite.set.size) {
      metabolite.set.tags[j] = trimws(met.line[j + 2])
    }
    existing.set = is.element(metabolite.set.tags, metabolite.labels)
    set.size = length(existing.set[existing.set == T])
    if ((set.size < met.size.threshold.min) || (set.size > met.size.threshold.max)) next
    temp.size.G[met.count] = set.size
    gs[met.count,] = c(metabolite.set.tags[existing.set], rep("null", max.size.G - temp.size.G[met.count]))
    temp.names[met.count] = metabolite.set.name
    temp.desc[met.count] = metabolite.set.desc
    met.count = met.count + 1
  }
  Ng = met.count - 1
  met.names = vector(length = Ng, mode = "character")
  met.desc = vector(length = Ng, mode = "character")
  size.G = vector(length = Ng, mode = "numeric")
  met.names = temp.names[1:Ng]
  met.desc = temp.desc[1:Ng]
  size.G = temp.size.G[1:Ng]
  
  N = length(A[,1])
  Ns = length(A[1,])
  print(c("Number of metabolites in dataset:", N))
  print(sprintf("Number of Metabolite Pathway Sets between size %d-%d: %d", met.size.threshold.min, met.size.threshold.max, Ng))
  print(c("Number of samples:", Ns))
  print(c("Original number of Metabolite Pathway Sets:", max.Ng))
  print(c("Maximum metabolite pathway set size:", max.size.G))
  
  # Read metabolite and metabolite set annotations if metabolite annotation file was provided
  all.metabolite.descs = vector(length = N, mode ="character")
  all.metabolite.symbols = vector(length = N, mode ="character")
  for (i in 1:N) {
    all.metabolite.descs[i] = metabolite.labels[i]
    all.metabolite.symbols[i] = metabolite.labels[i]
  }
  Obs.indicator = matrix(nrow= Ng, ncol=N)
  Obs.RES = matrix(nrow= Ng, ncol=N)
  Obs.ES = vector(length = Ng, mode = "numeric")
  Obs.arg.ES = vector(length = Ng, mode = "numeric")
  Obs.ES.norm = vector(length = Ng, mode = "numeric")
  
  # MSEA methodology
  # Compute observed and random permutation metabolite rankings
  obs.s2n = vector(length=N, mode="numeric")
  signal.strength = vector(length=Ng, mode="numeric")
  tag.frac = vector(length=Ng, mode="numeric")
  metabolite.frac = vector(length=Ng, mode="numeric")
  coherence.ratio = vector(length=Ng, mode="numeric")
  obs.phi.norm = matrix(nrow = Ng, ncol = nperm)
  correl.matrix = matrix(nrow = N, ncol = nperm)
  obs.correl.matrix = matrix(nrow = N, ncol = nperm)
  order.matrix = matrix(nrow = N, ncol = nperm)
  obs.order.matrix = matrix(nrow = N, ncol = nperm)
  
  nperm.per.call = 100
  n.groups = nperm %/% nperm.per.call
  n.rem = nperm %% nperm.per.call
  n.perms = c(rep(nperm.per.call, n.groups), n.rem)
  n.ends = cumsum(n.perms)
  n.starts = n.ends - n.perms + 1
  if (n.rem == 0) {n.tot = n.groups} else {n.tot = n.groups + 1}
  
  for (nk in 1:n.tot) {
    call.nperm = n.perms[nk]
    print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", 
                n.starts[nk], "--", n.ends[nk], sep=" "))
    O = MSEA.MetaboliteRanking(A, class.labels, metabolite.labels, call.nperm, sigma.correction = "MetaboliteCluster", replace=FALSE)
    order.matrix[,n.starts[nk]:n.ends[nk]] = O$order.matrix
    obs.order.matrix[,n.starts[nk]:n.ends[nk]] = O$obs.order.matrix
    correl.matrix[,n.starts[nk]:n.ends[nk]] = O$s2n.matrix
    obs.correl.matrix[,n.starts[nk]:n.ends[nk]] = O$obs.s2n.matrix
    rm(O)
  }
  obs.s2n = apply(obs.correl.matrix, 1, function(i) median(na.omit(i)))  # using median to assign enrichment scores
  obs.index = order(obs.s2n, decreasing=T)
  obs.s2n   = sort(obs.s2n, decreasing=T, na.last = TRUE)
  
  obs.metabolite.labels = metabolite.labels[obs.index]
  obs.metabolite.descs = all.metabolite.descs[obs.index]
  obs.metabolite.symbols = all.metabolite.symbols[obs.index]
  for (r in 1:nperm) {correl.matrix[, r] = correl.matrix[order.matrix[,r], r]}
  for (r in 1:nperm) {obs.correl.matrix[, r] = obs.correl.matrix[obs.order.matrix[,r], r]}
  metabolite.list2 = obs.index
  for (i in 1:Ng) {
    print(paste("Computing observed enrichment for metabolite set:", i, met.names[i], sep=" "))
    metabolite.set = gs[i,gs[i,] != "null"]
    metabolite.set2 = vector(length=length(metabolite.set), mode = "numeric")
    metabolite.set2 = match(metabolite.set, metabolite.labels)
    MSEA.results = MSEA.EnrichmentScore2(metabolite.list=metabolite.list2, metabolite.set=metabolite.set2, 
                                         weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
    Obs.ES[i] = MSEA.results$ES
    if (Obs.ES[i] >= 0) {  # compute signal strength
      tag.frac[i] = sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.G[i]
      metabolite.frac[i] = Obs.arg.ES[i]/N
    } else {
      tag.frac[i] = sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
      metabolite.frac[i] = (N - Obs.arg.ES[i] + 1)/N
    }
    signal.strength[i] = tag.frac[i] * (1 - metabolite.frac[i]) * (N / (N - size.G[i]))
  }
  
  # Compute enrichment for random permutations, using sample-label permutation
  phi = matrix(nrow = Ng, ncol = nperm)
  phi.norm = matrix(nrow = Ng, ncol = nperm)
  obs.phi = matrix(nrow = Ng, ncol = nperm)
  for (i in 1:Ng) {
    print(paste("Computing random permutations' enrichment for metabolite set:", i, met.names[i], sep=" "))
    metabolite.set = gs[i,gs[i,] != "null"]
    metabolite.set2 = vector(length=length(metabolite.set), mode = "numeric")
    metabolite.set2 = match(metabolite.set, metabolite.labels)
    for (r in 1:nperm) {
      metabolite.list2 = order.matrix[,r]
      # Fast version of enrichment scoring method
      MSEA.results = MSEA.EnrichmentScore2(metabolite.list=metabolite.list2, 
                                           metabolite.set=metabolite.set2, 
                                           weighted.score.type=weighted.score.type, 
                                           correl.vector=correl.matrix[, r])
      phi[i, r] = MSEA.results$ES
    }
    obs.metabolite.list2 = obs.order.matrix[,1]
    MSEA.results = MSEA.EnrichmentScore2(metabolite.list=obs.metabolite.list2, metabolite.set=metabolite.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
    obs.phi[i, 1] = MSEA.results$ES
    for (r in 2:nperm) {obs.phi[i, r] = obs.phi[i, 1]}
    gc()
  } 
  
  # Compute 3 types of p-values
  # Find nominal p-values
  print("Computing nominal p-values...")
  p.vals = matrix(1, nrow = Ng, ncol = 2)
  for (i in 1:Ng) {
    pos.phi = NULL
    neg.phi = NULL
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        pos.phi = c(pos.phi, phi[i, j])
      } else {
        neg.phi = c(neg.phi, phi[i, j])
      }
    }
    ES.value = Obs.ES[i]
    if (ES.value >= 0) {
      p.vals[i, 1] = signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=5)
    } else {
      p.vals[i, 1] = signif(sum(neg.phi <= ES.value)/length(neg.phi), digits=5)
    }
  }
  
  # Find effective size
  erf = function (x) {2 * pnorm(sqrt(2) * x)}
  KS.mean = function(N) { # KS mean as a function of set size N
    S = 0
    for (k in -100:100) {
      if (k == 0) next
      S = S + 4 * (-1)**(k + 1) * (0.25 * exp(-2 * k * k * N) - sqrt(2 * pi) *  erf(sqrt(2 * N) * k)/(16 * k * sqrt(N)))
    }
    return(abs(S))
  }
  # Rescaling normalization for each metabolite set null
  print("Computing rescaling normalization for each metabolite set null...")
  for (i in 1:Ng) {
    pos.phi = NULL
    neg.phi = NULL
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        pos.phi = c(pos.phi, phi[i, j])
      } else {
        neg.phi = c(neg.phi, phi[i, j])
      }
    }
    pos.m = mean(pos.phi)
    neg.m = mean(abs(as.numeric(neg.phi)))
    
    pos.phi = pos.phi/pos.m
    neg.phi = neg.phi/neg.m
    for (j in 1:nperm) {
      if (phi[i, j] >= 0) {
        phi.norm[i, j] = phi[i, j]/pos.m
      } else {
        phi.norm[i, j] = phi[i, j]/neg.m
      }
    }
    for (j in 1:nperm) {
      if (obs.phi[i, j] >= 0) {
        obs.phi.norm[i, j] = obs.phi[i, j]/pos.m
      } else {
        obs.phi.norm[i, j] = obs.phi[i, j]/neg.m
      }
    }
    if (Obs.ES[i] >= 0) {
      Obs.ES.norm[i] = Obs.ES[i]/pos.m
    } else {
      Obs.ES.norm[i] = Obs.ES[i]/neg.m
    }
  }
  
  # Compute FWER p-vals
  print("Computing FWER p-values...")
  max.ES.vals.p = NULL
  max.ES.vals.n = NULL
  for (j in 1:nperm) {
    pos.phi = NULL
    neg.phi = NULL
    for (i in 1:Ng) {
      if (phi.norm[i, j] >= 0) {
        pos.phi = c(pos.phi, phi.norm[i, j])
      } else {
        neg.phi = c(neg.phi, phi.norm[i, j])
      }
    }
    if (length(pos.phi) > 0) {
      max.ES.vals.p = c(max.ES.vals.p, max(pos.phi))
    }
    if (length(neg.phi) > 0) {
      max.ES.vals.n = c(max.ES.vals.n, min(neg.phi))
    }
  }
  for (i in 1:Ng) {
    ES.value = Obs.ES.norm[i]
    if (Obs.ES.norm[i] >= 0) {
      p.vals[i, 2] = signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), digits=5)
    } else {
      p.vals[i, 2] = signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), digits=5)
    }
  }
  
  # Compute FDRs
  print("Computing FDR q-values...")
  NES = vector(length=Ng, mode="numeric")
  phi.norm.mean  = vector(length=Ng, mode="numeric")
  obs.phi.norm.mean  = vector(length=Ng, mode="numeric")
  phi.norm.mean  = vector(length=Ng, mode="numeric")
  obs.phi.mean  = vector(length=Ng, mode="numeric")
  FDR.mean = vector(length=Ng, mode="numeric")
  
  Obs.ES.index = order(Obs.ES.norm, decreasing=T)
  Orig.index = seq(1, Ng)
  Orig.index = Orig.index[Obs.ES.index]
  Orig.index = order(Orig.index, decreasing=F)
  Obs.ES.norm.sorted = Obs.ES.norm[Obs.ES.index]
  met.names.sorted = met.names[Obs.ES.index]
  
  for (k in 1:Ng) {
    NES[k] = Obs.ES.norm.sorted[k]
    ES.value = NES[k]
    count.col = vector(length=nperm, mode="numeric")
    obs.count.col = vector(length=nperm, mode="numeric")
    for (i in 1:nperm) {
      phi.vec = phi.norm[,i]
      obs.phi.vec = obs.phi.norm[,i]
      if (ES.value >= 0) {
        count.col.norm = sum(phi.vec >= 0)
        obs.count.col.norm = sum(obs.phi.vec >= 0)
        count.col[i] = ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
        obs.count.col[i] = ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
      } else {
        count.col.norm = sum(phi.vec < 0)
        obs.count.col.norm = sum(obs.phi.vec < 0)
        count.col[i] = ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
        obs.count.col[i] = ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
      }
    }
    phi.norm.mean[k] = mean(count.col)
    obs.phi.norm.mean[k] = mean(obs.count.col)
    FDR.mean[k] = ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
  }
  FDR.mean[which(is.na(FDR.mean))] = 1

  # adjust q-values
  if (adjust.FDR.q.val == T) {
    pos.nes = length(NES[NES >= 0])
    min.FDR.mean = FDR.mean[pos.nes]
    for (k in seq(pos.nes - 1, 1, -1)) {
      if (FDR.mean[k] < min.FDR.mean) {min.FDR.mean = FDR.mean[k]}
      if (min.FDR.mean < FDR.mean[k]) {FDR.mean[k] = min.FDR.mean}
    }
    neg.nes = pos.nes + 1
    min.FDR.mean = FDR.mean[neg.nes]
    for (k in seq(neg.nes + 1, Ng)) {
      if (FDR.mean[k] < min.FDR.mean) {min.FDR.mean = FDR.mean[k]}
      if (min.FDR.mean < FDR.mean[k]) {FDR.mean[k] = min.FDR.mean}
    }
  }
  obs.phi.norm.mean.sorted = obs.phi.norm.mean[Orig.index]
  phi.norm.mean.sorted = phi.norm.mean[Orig.index]
  FDR.mean.sorted = FDR.mean[Orig.index]

  # Compute global statistic
  glob.p.vals = vector(length=Ng, mode="numeric")
  NULL.pass = vector(length=nperm, mode="numeric")
  OBS.pass = vector(length=nperm, mode="numeric")
  
  for (k in 1:Ng) {
    NES[k] = Obs.ES.norm.sorted[k]
    if (NES[k] >= 0) {
      for (i in 1:nperm) {
        NULL.pos = sum(phi.norm[,i] >= 0)
        NULL.pass[i] = ifelse(NULL.pos > 0, sum(phi.norm[,i] >= NES[k])/NULL.pos, 0)
        OBS.pos = sum(obs.phi.norm[,i] >= 0)
        OBS.pass[i] = ifelse(OBS.pos > 0, sum(obs.phi.norm[,i] >= NES[k])/OBS.pos, 0)
      }
    } else {
      for (i in 1:nperm) {
        NULL.neg = sum(phi.norm[,i] < 0)
        NULL.pass[i] = ifelse(NULL.neg > 0, sum(phi.norm[,i] <= NES[k])/NULL.neg, 0)
        OBS.neg = sum(obs.phi.norm[,i] < 0)
        OBS.pass[i] = ifelse(OBS.neg > 0, sum(obs.phi.norm[,i] <= NES[k])/OBS.neg, 0)
      }
    }
    glob.p.vals[k] = sum(NULL.pass >= mean(OBS.pass))/nperm
  }
  glob.p.vals.sorted = glob.p.vals[Orig.index]
  
  # Produce results report
  Obs.ES.norm = signif(Obs.ES.norm, digits=4)
  p.vals = signif(p.vals, digits=3)
  FDR.mean.sorted = signif(FDR.mean.sorted, digits=3)
  glob.p.vals.sorted = signif(glob.p.vals.sorted, digits=3)
  met.names = as.character(sapply(met.names, function(i) gsub("\\[1\\] ", "", i)))
  
  report = data.frame(Pathway=met.names, Size=size.G, NES=Obs.ES.norm, NOM.pval=p.vals[,1], FDR.qval=FDR.mean.sorted, FWER.pval=p.vals[,2], 
                      Global.pval=glob.p.vals.sorted)
  report = report[order(report$NES, decreasing=T),]
  report = report[which(report$NOM.pval<0.25),]
  colnames(report) = c("Pathway", "Size", "NES", "NOM\npval", "FDR\nqval", "FWER\npval", "Global\npval")
  
  return(report)
}

