# point transect testing

library(mrds)
library(testthat)

par.tol <- 1e-6
lnl.tol <- 1e-4

context("Point transects: ptexample")

test_that("Point transect example from distance gives same results",{

  data(ptdata.distance)
  xx <- ddf(dsmodel=~cds(key="hn", formula=~1), data=ptdata.distance,
            method="ds",
            meta.data=list(point=TRUE),truncation=max(ptdata.distance$distance))

  expect_that(xx$par,equals(2.283007,tol=par.tol))
  expect_that(xx$lnl,equals(-458.5701,tol=lnl.tol))
  expect_that(summary(xx)$average.p,equals(0.1644288,tol=par.tol))


})
