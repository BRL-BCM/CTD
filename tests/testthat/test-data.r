# stat.entropyFunction.r
test_that("Entropy function works", {
  expect_equal(stat.entropyFunction(c(1,1,1,0,0,0)), 1)
  expect_equal(stat.entropyFunction(c(0,0,0,0,0,0)), 0)
  expect_equal(stat.entropyFunction(c(1,1,1,1,1,1)), 0)
})

# stat.fishersMethod.r
test_that("Fisher's method works", {
  expect_equal(stat.fishersMethod(c(0, 0, 0)), 0)
  expect_equal(stat.fishersMethod(c(1, 1, 1)), 1)
  expect_equal(is.numeric(stat.fishersMethod(c(0.24, 0.34, 0.42))), TRUE)
})

# data.cohorts_coded.r
# data.Miller2015.r
# data.Wangler2017.r
# data.Thistlethwaite2020.r

# Z-score data
test_that("Zscore data works", {
  # data.zscoreData.r
  dis_data = matrix(rexp(500), ncol=100)
  rownames(dis_data) = sprintf("Feature%d", seq_len(nrow(dis_data)))
  colnames(dis_data) = sprintf("Sample%d", seq_len(ncol(dis_data)))
  ref_data = matrix(rexp(500), ncol=100)
  rownames(ref_data) = sprintf("Feature%d", seq_len(nrow(ref_data)))
  colnames(ref_data) = sprintf("Sample%d", seq_len(ncol(ref_data)))
  zscored.data = data.zscoreData(dis_data, ref_data)
  # Test1. surr.data should be a numeric matrix.
  expect_equal(class(zscored.data)[1], "matrix")
  expect_equal(all(apply(zscored.data, c(1,2), is.numeric)), TRUE)
  # Test2. If data.imputeData.r works, there should be no NAs.  
  expect_equal(any(is.na(zscored.data)), FALSE)
})

# Combine datasets.
test_that("Combine data works", {
  # data.combineData.r
  # Row names and column names are required for both input matrices.
  curr_data = matrix(rnorm(500), ncol=100)
  rownames(curr_data) = sprintf("Feature%d", sample(seq_len(20), 
                                              nrow(curr_data), 
                                            replace = FALSE))
  colnames(curr_data) = sprintf("Curr.Sample%d", seq_len(ncol(curr_data)))
  more_data = matrix(rnorm(500), ncol=100)
  rownames(more_data) = sprintf("Feature%d", sample(seq_len(20), 
                                            nrow(curr_data), 
                                            replace = FALSE))
  colnames(more_data) = sprintf("More.Sample%d", seq_len(ncol(curr_data)))
  combined.data = data.combineData(curr_data, more_data)
  expect_equal(class(combined.data)[1], "matrix")
  # Test1. Number of returned columns should be larger than or 
  #        equal to max number of columns of between the two input matrices.
  expect_equal(ncol(combined.data) >= max(ncol(curr_data), ncol(more_data)), TRUE)
  # Test2. Number of returned columns should be less than or 
  #        equal to sum of columns of between the two input matrices.
  expect_equal(ncol(combined.data) <= ncol(curr_data)+ncol(more_data), TRUE)
})

# Surrogate profile generation.
test_that("Surrogate profiles works", {
  curr_data = matrix(rnorm(1000), ncol=20)
  rownames(curr_data) = sprintf("Feature%d", seq_len(nrow(curr_data)))
  colnames(curr_data) = sprintf("Curr.Sample%d", seq_len(ncol(curr_data)))
  more_data = matrix(rnorm(1000), ncol=20)
  rownames(more_data) = sprintf("Feature%d", seq_len(nrow(more_data)))
  colnames(more_data) = sprintf("More.Sample%d", seq_len(ncol(more_data)))
  # Add missingness so data.imputeData.r will be tested.
  curr_data[sample(1:nrow(curr_data), 25, replace=FALSE), 
            sample(1:ncol(curr_data), 10, replace=FALSE)] = NA
  more_data[sample(1:nrow(more_data), 25, replace=FALSE), 
            sample(1:ncol(more_data), 10, replace=FALSE)] = NA
  surr.data = data.surrogateProfiles(curr_data, 1, ref_data = more_data)
  # Test1. surr.data should be a numeric matrix.
  expect_equal(class(surr.data)[1], "matrix")
  expect_equal(all(apply(surr.data, c(1,2), is.numeric)), TRUE)
  # Test2. If data.imputeData.r works, there should be no NAs.  
  expect_equal(any(is.na(surr.data)), FALSE)
  # Test3. If ref_data is specified, you should have both "disease_surr"
  #        and "control_surr" in the colnames of the returned matrix.
  expect_equal(any(grepl("disease_surr", colnames(surr.data))) && 
               any(grepl("control_surr", colnames(surr.data))), TRUE)
})

# Network pruning
test_that("Network pruning works", {
  ig = graph.adjacency(get.adjacency(sample_gnp(10, 0.8)), mode = "undirected", weighted = TRUE)
  ig_ref = graph.adjacency(get.adjacency(sample_gnp(10, 0.8)), mode = "undirected", weighted = TRUE)
  V(ig)$name = LETTERS[1:10]
  V(ig_ref)$name = LETTERS[1:10]
  # graph.naivePruning.r
  ig_pruned = graph.naivePruning(ig, ig_ref)
  expect_equal(length(E(ig_pruned)$weight) <= length(E(ig)$weight), TRUE)
})

# Probability diffusion algorithm
test_that("Probability diffusion works", {
  adj_mat = rbind(c(0,3,1,0,0,0,0,0,0), #A's neighbors
                  c(3,0,2,2,0,0,0,0,0), #B's neighbors
                  c(1,2,0,0,2,0,0,0,0), #C's neighbors
                  c(0,2,0,0,1,0,1,1,0), #D's neighbors
                  c(0,0,2,1,0,2,0,2,0), #E's neighbors
                  c(0,0,0,0,2,0,0,0,0), #F's neighbors
                  c(0,0,0,1,0,0,0,1,0), #G's neighbors
                  c(0,0,0,1,2,0,1,0,1), #H's neighbors
                  c(0,0,0,0,0,0,0,1,0) #I's neighbors
  )
  rownames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  colnames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  G = vector("numeric", length=ncol(adj_mat))
  names(G)=colnames(adj_mat)
  ig = graph.adjacency(adj_mat, mode="undirected", weighted=TRUE, add.colnames = "name")
  coords = layout.fruchterman.reingold(ig)
  .GlobalEnv$imgNum = 1
  G_new = graph.diffuseP1(p1=1.0, sn="A", G=G, vNodes="A", 
                          thresholdDiff=0.01, adj_mat=adj_mat, verbose=FALSE)
  # Inherited probabilities across all nodes should add to 1.
  # Which node inherited the highest probability from startNode?
  expect_equal(names(G_new[which.max(G_new)]), "B")
  # Test1. Returned probability associated with start node should be 0.
  expect_equal(G_new[["A"]], 0)
  # Test2. Sum of returned probabilities over all nodes should be 1.
  expect_equal(sum(unlist(G_new)), 1)
  # Test3. No errors thrown when movie is generated (means graph.takeDiffusionSnapShot.r works)
  expect_equal(length(G_new)==length(G), TRUE)
  # Test4. No errors thrown when start node is stranded by its visited nodes 
  #        (means graph.connectToExt.r works).
  G_new = graph.diffuseP1(p1=1.0, sn="A", G=G, vNodes=c("A","B","C"), 
                          thresholdDiff=0.01, adj_mat=adj_mat)
  expect_equal(length(G_new)==length(G), TRUE)
})

# Multi-node encoding
test_that("Multi-node pipeline works", {
  adj_mat = rbind(c(0,3,1,0,0,0,0,0,0), #A's neighbors
                  c(3,0,2,2,0,0,0,0,0), #B's neighbors
                  c(1,2,0,0,2,0,0,0,0), #C's neighbors
                  c(0,2,0,0,1,0,1,1,0), #D's neighbors
                  c(0,0,2,1,0,2,0,2,0), #E's neighbors
                  c(0,0,0,0,2,0,0,0,0), #F's neighbors
                  c(0,0,0,1,0,0,0,1,0), #G's neighbors
                  c(0,0,0,1,2,0,1,0,1), #H's neighbors
                  c(0,0,0,0,0,0,0,1,0) #I's neighbors
  )
  rownames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  colnames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  G = vector("numeric", length=ncol(adj_mat))
  names(G)=colnames(adj_mat)
  # multiNode.getNodeRanks.r, with movie functionality
  ig = graph.adjacency(adj_mat, mode="undirected", weighted=TRUE, add.colnames = "name")
  coords = layout.fruchterman.reingold(ig)
  S = names(G)[1:3]
  ranks = multiNode.getNodeRanks(S=S, G=G, p1=0.9, thresholdDiff=0.01, adj_mat=adj_mat,
                                 num.misses=log2(length(G)), verbose=FALSE)
  # If no error, we know graph.takeNetWalkSnapShot.r is also working.
  expect_equal(all(lapply(ranks, length)>0), TRUE)
  expect_equal(all(unlist(lapply(ranks, function(i) any(S %in% i)))), TRUE)
  # mle.getPtBSbyK.r
  ptBSbyK = mle.getPtBSbyK(S, ranks)
  expect_equal(all(lapply(ptBSbyK, sum)>0), TRUE)
  # mle.getEncodingLength.r
  res = mle.getEncodingLength(ptBSbyK, NULL, NULL, G)
  expect_equal(nrow(res)==length(S), TRUE)
  expect_equal(max(res$d.score)>0, TRUE)
})

# Single-node encoding
test_that("Single node pipeline works", {
  adj_mat = rbind(c(0,3,1,0,0,0,0,0,0), #A's neighbors
                  c(3,0,2,2,0,0,0,0,0), #B's neighbors
                  c(1,2,0,0,2,0,0,0,0), #C's neighbors
                  c(0,2,0,0,1,0,1,1,0), #D's neighbors
                  c(0,0,2,1,0,2,0,2,0), #E's neighbors
                  c(0,0,0,0,2,0,0,0,0), #F's neighbors
                  c(0,0,0,1,0,0,0,1,0), #G's neighbors
                  c(0,0,0,1,2,0,1,0,1), #H's neighbors
                  c(0,0,0,0,0,0,0,1,0) #I's neighbors
  )
  rownames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  colnames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  G = vector("numeric", length=ncol(adj_mat))
  names(G)=colnames(adj_mat)
  # singleNode.getNodeRanksN.r, with movie functionality
  ig = graph.adjacency(adj_mat, mode="undirected", weighted=TRUE, add.colnames = "name")
  coords = layout.fruchterman.reingold(ig)
  S = names(G)[1:3]
  ranks=list()
  # If S is not specified, ranks should be length of G
  for (i in 1:length(S)) {
    ranks[[i]] = singleNode.getNodeRanksN(n=i, G=G, p1=0.9, thresholdDiff=0.01, adj_mat=adj_mat, 
                                          S=S, num.misses=log2(length(G)), verbose=FALSE)
  }
  names(ranks)=S
  # If no error, we know graph.takeNetWalkSnapShot.r is also working.
  expect_equal(all(lapply(ranks, length)>0), TRUE)
  # Try without specifying S. Should be length of G
  ranks[[1]] = singleNode.getNodeRanksN(n=1, G=G, p1=0.9, thresholdDiff=0.01, adj_mat=adj_mat)
  expect_equal(length(ranks[[1]])==length(G), TRUE)
  # mle.getPtBSbyK.r
  ptBSbyK = mle.getPtBSbyK(S, ranks)
  expect_equal(all(lapply(ptBSbyK, sum)>0), TRUE)
  # mle.getEncodingLength.r
  res = mle.getEncodingLength(ptBSbyK, NULL, NULL, G)
  expect_equal(nrow(res)==length(S), TRUE)
  expect_equal(max(res$d.score)>0, TRUE)
})

# Patient distances 
test_that("Patient distances works", {
  adj_mat = rbind(c(0,3,1,0,0,0,0,0,0), #A's neighbors
                  c(3,0,2,2,0,0,0,0,0), #B's neighbors
                  c(1,2,0,0,2,0,0,0,0), #C's neighbors
                  c(0,2,0,0,1,0,1,1,0), #D's neighbors
                  c(0,0,2,1,0,2,0,2,0), #E's neighbors
                  c(0,0,0,0,2,0,0,0,0), #F's neighbors
                  c(0,0,0,1,0,0,0,1,0), #G's neighbors
                  c(0,0,0,1,2,0,1,0,1), #H's neighbors
                  c(0,0,0,0,0,0,0,1,0) #I's neighbors
  )
  rownames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  colnames(adj_mat) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  data_mx = apply(Miller2015[-1,which(Miller2015[1,]=="Argininemia")], c(1,2), as.numeric)
  data_mx = data_mx[sample(1:nrow(data_mx), nrow(adj_mat), replace=FALSE),]
  rownames(data_mx) = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
  G = vector("numeric", length=ncol(adj_mat))
  names(G)=colnames(adj_mat)
  ranks = list()
  for (n in seq_len(length(G))) { 
   ranks[[n]] = singleNode.getNodeRanksN(n, G, p1=0.9, thresholdDiff=0.01, adj_mat) 
  }
  names(ranks) = names(G)
  # Also pre-compute patient bitstrings for faster distance calculations.
  ptBSbyK = list()
  for (pt in seq_len(ncol(data_mx))) {
   S = names(sort(abs(data_mx[,pt]),decreasing=TRUE)[seq_len(3)])
   ptBSbyK[[colnames(data_mx)[pt]]] = mle.getPtBSbyK(S, ranks)
  }
  # Build your results ("res") list object to store patient distances at different size k's.
  res = list()
  t = list(ncd=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
  colnames(t$ncd) = colnames(data_mx)
  rownames(t$ncd) = colnames(data_mx)
  for (i in seq_len(3)) { res[[i]] = t }
  for (pt in seq_len(ncol(data_mx))) {
    ptID = colnames(data_mx)[pt]
    for (pt2 in pt:ncol(data_mx)) {
      ptID2 = colnames(data_mx)[pt2]
      # mle.getPtDist.r
      tmp = mle.getPtDist(ptBSbyK[[ptID]], ptID, ptBSbyK[[ptID2]], ptID2, data_mx, ranks, p1=0.9, thresholdDiff=0.01, adj_mat)
      for (k in seq_len(3)) {
        res[[k]]$ncd[ptID, ptID2] = tmp$NCD[k]
        res[[k]]$ncd[ptID2, ptID] = tmp$NCD[k]
      }
    }
  }
  expect_equal(length(res)==3, TRUE)
  expect_equal(unlist(lapply(res, function(i) any(is.na(i$ncd)))), c(FALSE, FALSE, FALSE))
  expect_equal(all(unlist(lapply(res, function(i) all(i$ncd>=0)))), TRUE)
  expect_equal(all(unlist(lapply(res, function(i) all(diag(i$ncd)==0)))), TRUE)
  # mle.getMinPtDistance.r
  min_dist = mle.getMinPtDistance(lapply(res, function(i) i$ncd))
  expect_equal(all(min_dist>=0) && all(min_dist<=1), TRUE)
})

















