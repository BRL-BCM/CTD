# stat.entropyFunction.r
test_that("Entropy function works", {
  expect_equal(stat.entropyFunction(c(1,1,1,0,0,0)), 1)
  expect_equal(stat.entropyFunction(c(0,0,0,0,0,0)), 0)
  expect_equal(stat.entropyFunction(c(1,1,1,1,1,1)), 0)
})

# stat.fishersMethod.r
test_that("Fisher's method works", {
  expect_equal(stat.fishersMethod(c(0.05, 0.25, 1.0)), 0.1872888)
})

# data.cohorts_coded.r
# data.Miller2015.r
# data.Wangler2017.r
# data.Thistlethwaite2020.r


# data.imputeData.r

# data.combineData.r

# data.surrogateProfiles.r

# data.zscoreData.r

# graph.connectToExt.r

# graph.diffuseP1.r

# graph.naivePruning.r

# graph.takeDiffusionSnapShot.r

# graph.takeNetWalkSnapShot.r

# mle.getEncodingLength.r

# mle.getMinPtDistance.r

# mle.getPtDist.r

# multiNode.getNodeRanks.r

# multiNode.getPtBSbyK.r

# singleNode.getNodeRanksN.r

# singleNode.getPtBSbyK.r









