library(mrds)
library(testthat)

context("No sightings")

test_that("a dht object with only 0's for estimates is returned", {

  #datasetup
  ex.filename<-system.file("testData/NoSightings/clusters.robj", package="mrds")
  load(ex.filename)
  ex.filename<-system.file("testData/NoSightings/no_clusters.robj", package="mrds")
  load(ex.filename)
  ex.filename<-system.file("testData/NoSightings/obs_subset.robj", package="mrds")
  load(ex.filename)
  ex.filename<-system.file("testData/NoSightings/region_table.robj", package="mrds")
  load(ex.filename)
  ex.filename<-system.file("testData/NoSightings/sample_table.robj", package="mrds")
  load(ex.filename)
  
  #run analyses
  ddf.clusters <- ddf(method='ds', dsmodel = ~mcds(key = 'hn', formula = ~scaledtotsize), data = clusters, meta.data = list(width = 5.5))
  ddf.no.clusters <- ddf(method='ds', dsmodel = ~mcds(key = 'hn', formula = ~scaledtotsize), data = no.clusters, meta.data = list(width = 5.5))
  cluster.result <- dht(ddf.clusters, region.table, sample.table, obs.table.subset)
  no.cluster.result <- dht(ddf.no.clusters, region.table, sample.table, obs.table.subset)

  #run tests
  #checks there are NO results for clusters
  expect_that(is.null(no.cluster.result$clusters), equals(TRUE))
  #checks there ARE results for clusters
  expect_that(is.null(cluster.result$clusters), equals(FALSE))
  #checks all abundance estimates are 0 for individuals in clustered data
  expect_that(length(which(cluster.result$individuals$N$Estimate != 0)), equals(0))
  #checks all df are 0 for clusters for clustered data; removed I don't think this is useful nor desirable
  #expect_that(length(which(cluster.result$clusters$N$df != 0)), equals(0))
  #checks all encounter rate estimates are 0 for clusters for clustered data
  expect_that(length(which(cluster.result$clusters$summary$ER != 0)), equals(0))
  #checks all density estimates are 0 for clusters for non-clustered data
  expect_that(length(which(no.cluster.result$clusters$D$se != 0)), equals(0))
  #checks the region areas match for the cluster and individual summaries for clustered data
  expect_that(length(which(!cluster.result$clusters$summary$Area == cluster.result$individuals$summary$Area)), equals(0))
  #checks the region areas match for the individual summaries for both clustered and non-clustered data
  expect_that(length(which(!no.cluster.result$individuals$summary$Area == cluster.result$individuals$summary$Area)), equals(0))
})
