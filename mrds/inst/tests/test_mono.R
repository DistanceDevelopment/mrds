# Monotonicity constraints test file
# test that using Rsolnp gives the same results
# as optimx() when it should!

library(testthat)
library(mrds)

# boring bookexamples-based testing...

context("Monotonicity: bookexamples tests")

test_that("parameter estimates are the same", { 

   # want to make sure that we get the same answers when
   # the function is monotone by default -- ie half-normal

   # data setup  
   data(book.tee.data)
   region<-book.tee.data$book.tee.region
   egdata<-book.tee.data$book.tee.dataframe
   samples<-book.tee.data$book.tee.samples
   obs<-book.tee.data$book.tee.obs

   # run using optimx
   result.optimx<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1), 
                      data = egdata[egdata$observer ==1, ], method = "ds", 
                      meta.data = list(width = 4))
   # run with Rsolnp, mono=TRUE, mono.strict=FALSE
   result.mono<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1), 
                    data = egdata[egdata$observer ==1, ], method = "ds", 
                    meta.data = list(width = 4,mono=TRUE,mono.strict=FALSE))
   # run with Rsolnp, mono=TRUE, mono.strict=TRUE
   result.strict<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1), 
                      data = egdata[egdata$observer ==1, ], method = "ds", 
                      meta.data = list(width = 4,mono=TRUE,mono.strict=TRUE))
#     summary(result,se=TRUE)


   # test that the parameter is the same
   expect_that(result.optimx$par, equals(result.mono$par,tolerance=1e-5)) # FAIL 
   expect_that(result.optimx$par, equals(result.strict$par,tolerance=1e-5)) # FAIL
   expect_that(result.mono$par, equals(result.strict$par,tolerance=1e-5))

   # check that the summary gives the right answers too... 
   # these won't work anymore since if we don't have adjustments, we just
   # switch to optimx()
   #expect_that(summary(result.mono),prints_text("Monotonicity constraints were enforced."))
   #expect_that(summary(result.strict),prints_text("Strict monotonicity constraints were enforced."))
   # check that the optimx() summary doesn't print anything about monotonicity
   expect_that(summary(result.optimx),prints_text("[^(Monotonicity constraints were enforced.)]"))
   expect_that(summary(result.optimx),prints_text("[^(Strict monotonicity constraints were enforced.)]"))
})

test_that("monotonicity warnings are correct", {
   # data setup  
   data(book.tee.data)
   region<-book.tee.data$book.tee.region
   egdata<-book.tee.data$book.tee.dataframe
   samples<-book.tee.data$book.tee.samples
   obs<-book.tee.data$book.tee.obs

   # if the user asks for monotonicity in a covariate model
   # we should warn and then use optimx()


   # check that there is a warning
   expect_that(result<-ddf(dsmodel = ~mcds(key = "hn", formula = ~sex), 
                   data = egdata[egdata$observer ==1, ], method = "ds", 
                   meta.data = list(width = 4,mono=TRUE)),
               gives_warning("Covariate models cannot be constrained for monotonicity."))

   # check that the model that was fit was without monotonicity
   expect_that(result$ds$aux$mono,equals(FALSE))
   expect_that(result$ds$aux$mono.strict,equals(FALSE))
   #check that the summary doesn't say anything about monotonicity
   expect_that(print(summary(result)),
               prints_text("[^(Monotonicity constraints were enforced.)]"))
   expect_that(print(summary(result)),
             prints_text("[^(Strict monotonicity constraints were enforced.)]"))
})





context("Monotonicity: mixture tests")
test_that("monotonic and non-monotonic fits are different for non-monotone data",{

   # save the seed
   #rngsave<-.Random.seed
   set.seed(3141)

   # simulate some non-monotonic data
   dat<-mrds:::sim.mix(100,c(0.1,3),c(0.3,0.7),10,means=c(0,4))
   dat<-data.frame(distance=dat,object=1:length(dat),observed=rep(1,length(dat)))

   trunc<-7

   # fit without constraint
   result.n<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2)), 
                 data=dat, method = "ds", 
                 meta.data=list(width=trunc,mono=FALSE))
   # with weak monotonicity
   result.w<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2)), 
                 data=dat, method = "ds", 
                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
   # with strong monotonicity
   result.s<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2)), 
                 data=dat, method = "ds", 
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


