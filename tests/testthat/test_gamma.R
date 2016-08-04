library(mrds)
library(testthat)

par.tol<-1e-3

context("gamma detection function")

test_that("standardisation is done correctly",{

  # generate some gamma data where the mode is basically 0
  set.seed(123)
  gdata<-data.frame(object=1:1000,distance=rgamma(1000,scale=1,shape=1))

  # this will complain about first partial hessian being singular
  # suppress that
  mm <- suppressWarnings(ddf(dsmodel=~cds(key="gamma"), data=gdata, method="ds",
                             meta.data=list(width=max(gdata$distance))))

  # check fitted parameters
  parcomp <- mm$par
  names(parcomp) <- NULL
  expect_equal(c(-19.097429433, -0.001956123),parcomp,tol=par.tol)


  # now with adjustments
  # again this will complain about first partial hessian being singular
  # suppress that
  mm <- suppressWarnings(ddf(dsmodel=~cds(key="gamma",adj.series="cos",
                                          adj.order=2),
                             data=gdata, method="ds",
                             meta.data=list(width=max(gdata$distance))))

  # check fitted parameters
  parcomp <- mm$par
  names(parcomp) <- NULL
  expect_equal(c(-21.6765623,0.1103963, 0.2100172),parcomp,tol=par.tol)

  # generate some gamma data where the mode is away from 0
  set.seed(123)
  gdata<-data.frame(object=1:1000,distance=rgamma(1000,scale=2,shape=3))

  mm <- ddf(dsmodel=~cds(key="gamma"), data=gdata, method="ds",
            meta.data=list(width=max(gdata$distance)))

  # check fitted parameters
  parcomp <- mm$par
  names(parcomp) <- NULL
  expect_equal(c(0.7517943,1.9497688),parcomp,tol=par.tol)

})
