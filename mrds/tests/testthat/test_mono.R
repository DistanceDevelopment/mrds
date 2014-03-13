# Monotonicity constraints test file
# test that using Rsolnp gives the same results
# as optimx() when it should!

library(testthat)
library(mrds)

tol <- 1e-6

# boring bookexamples-based testing...
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe

context("Monotonicity: bookexamples tests")

test_that("bookexamples parameter estimates and likelihood is correct", {

  # run with Rsolnp, mono=TRUE, mono.strict=FALSE
  result.mono<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1, adj.series="cos",
                                   adj.order=2),
                   data = egdata[egdata$observer ==1, ], method = "ds",
                   meta.data = list(width = 4,mono=TRUE,mono.strict=FALSE))

  mono.par <- c(0.663253549669784, -3.25134503272446e-07)
  names(result.mono$par) <- NULL
  expect_that(result.mono$par,equals(mono.par,tol=tol))
  expect_that(result.mono$lnl,equals(-154.569227292222,tol=tol))

  # run with Rsolnp, mono=TRUE, mono.strict=TRUE
  result.strict<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1, adj.series="cos",
                                     adj.order=2),
                     data = egdata[egdata$observer ==1, ], method = "ds",
                     meta.data = list(width = 4,mono=TRUE,mono.strict=TRUE))

  strict.par <- c(0.663257996011827,-3.27954455074367e-07)
  names(result.strict$par) <- NULL
  expect_that(result.strict$par,equals(strict.par,tol=tol))
  expect_that(result.strict$lnl,equals(-154.56922729565,tol=tol))

})

test_that("ddf when monotonicity is not required", {

  # covariate models die
  expect_that(result<-ddf(dsmodel = ~mcds(key = "hn", formula = ~sex),
                  data = egdata[egdata$observer ==1, ], method = "ds",
                  meta.data = list(width = 4,mono=TRUE)),
      throws_error("Covariate models cannot be constrained for monotonicity."))

  #key only models die
  expect_that(result<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1),
                  data = egdata[egdata$observer ==1, ], method = "ds",
                  meta.data = list(width = 4,mono=TRUE)),
     shows_message("Key only models do not require monotonicity contraints. Not constraining model for monotonicity."))

})


context("Monotonicity: mixture tests")

test_that("monotonic and non-monotonic fits differ for non-monotone data",{

   # save the seed
   #rngsave<-.Random.seed
   set.seed(3141)

   # simulate some non-monotonic data
   dat<-mrds:::sim.mix(100,c(0.1,3),c(0.3,0.7),10,means=c(0,4))
   dat<-data.frame(distance=dat,
                   object=1:length(dat),
                   observed=rep(1,length(dat)))

   trunc<-7

   # fit without constraint
   result.n<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",
                 adj.order=c(2)), data=dat, method = "ds",
                 meta.data=list(width=trunc,mono=FALSE))
   # with weak monotonicity
   result.w<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",
                 adj.order=c(2)), data=dat, method = "ds",
                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
   # with strong monotonicity
   result.s<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",
                 adj.order=c(2)), data=dat, method = "ds",
                 meta.data=list(width=trunc,mono=TRUE,mono.strict=TRUE))

   # plot
   #par(mfrow=c(1,3))
   #plot(result.n);plot(result.w);plot(result.s)

   # evaluate the detection function at a bunch of points
   x<-seq(0,trunc,len=1000)
   newdata<-data.frame(distance=x,object=1:length(x),observed=rep(1,length(x)))

   ddfobj<-result.n$ds$aux$ddfobj
   ddfobj$scale$dm<-mrds:::setcov(newdata, as.formula(ddfobj$scale$formula))$cov
   pred.n<-mrds:::detfct(x,ddfobj,width=trunc)

   ddfobj<-result.w$ds$aux$ddfobj
   ddfobj$scale$dm<-mrds:::setcov(newdata, as.formula(ddfobj$scale$formula))$cov
   pred.w<-mrds:::detfct(x,ddfobj,width=trunc)

   ddfobj<-result.s$ds$aux$ddfobj
   ddfobj$scale$dm<-mrds:::setcov(newdata, as.formula(ddfobj$scale$formula))$cov
   pred.s<-mrds:::detfct(x,ddfobj,width=trunc)

   expect_that(all(diff(order(pred.n,decreasing=TRUE))==1),is_false())
   expect_that(all(diff(order(pred.w,decreasing=TRUE))==1),is_true())
   expect_that(all(diff(order(pred.s,decreasing=TRUE))==1),is_true())

   # recover the seed
   #.Random.seed<-rngsave
})


