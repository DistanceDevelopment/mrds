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


context("Point transects: variance")

test_that("P3 variance estimator", {

  # use the pt example data
  data(ptdata.distance)

  # generate some effort data
  obs <- unique(ptdata.distance[, c("object", "Sample.Label", "Region.Label")])
  sample <- unique(ptdata.distance[, c("Sample.Label", "Region.Label")])
  region <- data.frame(Region.Label=1, Area=1)
  ptdata.distance$size <- 1

  # fit a model
  xx <- ddf(dsmodel=~cds(key="hn", formula=~1), data=ptdata.distance,
            method="ds",
            meta.data=list(point=TRUE,width=max(ptdata.distance$distance)))

  # with equal effort
  sample$Effort <- 1
  d <- dht(xx, obs.table=obs, region.table=region, sample.table=sample)
  expect_equal(d$individuals$N$se, 0.001000225155)

  # make effort unequal
  set.seed(1234)
  sample$Effort <- sample(1:3, nrow(sample), replace=TRUE)
  d <- dht(xx, obs.table=obs, region.table=region, sample.table=sample)
  expect_equal(d$individuals$N$se, 0.0007365085886)


  # save code for later for testing...
  ## pt encounter-rate variance "P3" from Fewster et al (2009)
  #tk <- sample$Effort
  #K <- nrow(sample)
  #T <- sum(tk)

  ## merge things
  #mm <- merge(region, sample, by="Region.Label")
  #mm <- merge(mm, obs, by="Sample.Label")
  #mm <- merge(mm, ptdata.distance, by="object")
  #mm$N <- mm$size/predict(xx)[[1]]
  #Nck <-  aggregate(mm$N, by=list(mm$Sample.Label.x), sum)$x
  #Nc <- sum(Nck)

  ## these are the same, if we calculate Nhatc rather than Nhat
  ## i.e., ignore the (A/a)^2 multiplier out front
  #ervar <- #(1/(T*pi*2*max(ptdata.distance$distance)))^2 *
  #          1/(T*(K-1)) * sum(tk * (Nck/tk - Nc/T)^2)

  #ltervar <- #(1/(2*max(ptdata.distance$distance)*T))^2*
  #           K/(T^2*(K-1)) * sum(tk^2 * (Nck/tk - Nc/T)^2)

  #varnn <- mrds:::varn(sample$Effort, Nck, type="R2")
  #varnn3 <- mrds:::varn(tk, Nck, type="R3")

  ## things go bad with uneven effort...
  #tk <- sample(1:3, K, replace=TRUE)
  #T <- sum(tk)

  #ervar <- #(1/(T*pi*2*max(ptdata.distance$distance)))^2 *
  #          1/(T*(K-1)) * sum(tk * (Nck/tk - Nc/T)^2)

  #ltervar <- #(1/(2*max(ptdata.distance$distance)*T))^2*
  #           K/(T^2*(K-1)) * sum(tk^2 * (Nck/tk - Nc/T)^2)

  #varnn <- mrds:::varn(tk, Nck, type="R2")
  #varnn3 <- mrds:::varn(tk, Nck, type="R3")


})
