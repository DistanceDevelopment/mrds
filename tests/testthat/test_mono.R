# Monotonicity constraints test file
# test that using Rsolnp gives the same results
# as optimx() when it should!

par.tol <- 1e-6
tol <- 1e-4

# boring bookexamples-based testing...
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe

context("Monotonicity: bookexamples tests")

test_that("bookexamples par. ests. and likelihood are correct", {

  # run with Rsolnp, mono=TRUE, mono.strict=FALSE
  result.mono<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1, adj.series="cos",
                                   adj.order=2),
                   data = egdata[egdata$observer ==1, ], method = "ds",
                   meta.data = list(width = 4,mono=TRUE,mono.strict=FALSE))

  mono.par <- c(0.66068519,  -0.01592319)
  names(result.mono$par) <- NULL
  expect_equal(result.mono$par, mono.par, tolerance=par.tol)
  expect_equal(result.mono$lnl, -154.56193, tol=tol)

  # run with Rsolnp, mono=TRUE, mono.strict=TRUE
  result.strict<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1, adj.series="cos",
                                     adj.order=2),
                     data = egdata[egdata$observer ==1, ], method = "ds",
                     meta.data = list(width = 4,mono=TRUE,mono.strict=TRUE))

  strict.par <- c(0.6606852,-0.0159232)
  names(result.strict$par) <- NULL
  expect_equal(result.strict$par,strict.par,tol=tol)
  expect_equal(result.strict$lnl,-154.56193,tol=tol)

})

test_that("ddf tells you when monotonicity is not required", {

  ## covariate models die
  expect_that(result<-ddf(dsmodel = ~mcds(key = "hn", formula = ~sex),
                  data = egdata[egdata$observer ==1, ], method = "ds",
                  meta.data = list(width = 4,mono=TRUE)),
      throws_error("Covariate models cannot be constrained for monotonicity."))

  #key only models die
  expect_that(result<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
                  data = egdata[egdata$observer ==1, ], method = "ds",
                  meta.data = list(width = 4,mono=TRUE)),
     shows_message("Key only model: not constraining for monotonicity."))

})


context("Monotonicity: mixture tests")

# save the seed
#rngsave<-.Random.seed
set.seed(341)

# simulate some non-monotonic data
dat<-mrds:::sim.mix(1000,c(0.1,0.5),c(0.2,0.8),10,means=c(0,2.5))
dat<-data.frame(distance=dat,
                object=1:length(dat),
                observed=rep(1,length(dat)))

trunc<-5

# fit without constraint
expect_warning(result.n<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,
                             adj.series="cos", adj.order=c(2,3)),
                             data=dat, method = "ds",
                             meta.data=list(width=trunc,mono=FALSE)),
               "Detection function is not strictly monotonic!")
# with weak monotonicity
expect_warning(result.w<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,
                             adj.series="cos", adj.order=c(2,3)),
                             data=dat, method = "ds",
                             meta.data=list(width=trunc,mono=TRUE,
                                            mono.strict=FALSE)),
               "Detection function is not strictly monotonic!")
# with strong monotonicity
result.s<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",
              adj.order=c(2,3)), data=dat, method = "ds",
              meta.data=list(width=trunc,mono=TRUE,mono.strict=TRUE))

## plot
#par(mfrow=c(2,3))
#plot(result.n,main="no constraints")
#plot(result.w,main="weak constraints")
#plot(result.s,main="strict constraints")
#check.mono(result.n,plot=TRUE,n.pts=20)
#check.mono(result.w,plot=TRUE,n.pts=20)
#check.mono(result.s,plot=TRUE,n.pts=20)

test_that("non-monotonic for non-monotone data",{
  # expect a non-monotonic fit
  expect_warning(mono.chk <- check.mono(result.n,n.pts=20),"Detection function is not strictly monotonic!")
  expect_that(mono.chk,is_false())
})


test_that("weakly monotone for weakly monotone constraints",{
  # check the weak fit is weak
  expect_that(check.mono(result.w,strict=FALSE,n.pts=20),is_true())

  expect_warning(mono.chk <- check.mono(result.w,n.pts=20),"Detection function is not strictly monotonic!")
  expect_that(mono.chk,is_false())
})


test_that("strictly monotonic for strictly monotone constaints",{
  # expect that the strict fit is strict
  expect_that(check.mono(result.s,strict=FALSE,n.pts=20),is_true())
  expect_that(check.mono(result.s,strict=TRUE,n.pts=20),is_true())
})



#context("Monotonicity checks for covariates")
#
#test_that("Montonicity checks for covariate models",{
#  result.nm<-ddf(dsmodel = ~mcds(key = "hn", formula = ~size, adj.series="cos",
#                                 adj.order=2,adj.scale="scale"),
#                 data = egdata[egdata$observer ==1, ], method = "ds",
#                 meta.data = list(width = 4))
#mrds:::check.mono(result.nm,plot=TRUE)
#
#
#
#
#
#})

