# point transect testing

library(mrds)
library(testthat)

par.tol <- 1e-5
lnl.tol <- 1e-4

context("Point transects: ptexample")

test_that("Point transect example from distance gives same results",{

  data(ptdata.distance)
  xx <- ddf(dsmodel=~cds(key="hn", formula=~1), data=ptdata.distance,
            method="ds",
            meta.data=list(point=TRUE,width=max(ptdata.distance$distance)))#,
            #control=list(upperbounds=Inf,lowerbounds=-Inf,showit=3))
            #control=list(optimx.maxit=10000,showit=3))

  expect_that(xx$par,equals(2.283007,tol=par.tol))
  expect_that(xx$lnl,equals(-458.5701,tol=lnl.tol))
  expect_that(summary(xx)$average.p,equals(0.1644288,tol=par.tol))


})


context("Point transects: simulated data")

test_that("Nothing has changed?", {

  # thanks to Tiago Marques for this example
  set.seed(123)
  n <- 1000
  w <- 8000
  sigma <- 3000
  xs <- runif(n,-w,w)
  ys <- runif(n,-w,w)
  rs <- sqrt(xs^2+ys^2)
  gx <- function(x,sigma){exp(-(x/sigma)^2)}
  px <- gx(rs,sigma)
  pdet <- runif(n,0,1)
  det <- px>pdet
  data4D <- data.frame(distance=rs[det==1])
  data4D$object <- 1:nrow(data4D)
  data4D$detected <- rep(1, nrow(data4D))


  # half-normal
  xx <- ddf(dsmodel=~cds(key="hn", formula=~1), data=data4D,
            method="ds",
            meta.data=list(point=TRUE, width=8000))

  expect_that(xx$par, equals(7.6438325, tol=par.tol))
  expect_that(xx$lnl, equals(-1066.120945, tol=lnl.tol))
  expect_that(summary(xx)$average.p, equals(0.1361185, tol=par.tol))


  # hazard-rate
  xx <- ddf(dsmodel=~cds(key="hr", formula=~1), data=data4D,
            method="ds",
            meta.data=list(point=TRUE, width=8000))

  expect_that(unname(xx$par), equals(c(1.554381, 7.943638), tol=par.tol))
  expect_that(xx$lnl, equals(-1069.143006, tol=lnl.tol))
  expect_that(summary(xx)$average.p, equals(0.1862252, tol=par.tol))



})
