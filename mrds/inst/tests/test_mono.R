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
   region<<-book.tee.data$book.tee.region
   egdata<<-book.tee.data$book.tee.dataframe
   samples<<-book.tee.data$book.tee.samples
   obs<<-book.tee.data$book.tee.obs

   # run using optimx
   result.optimx<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1), 
                      data = egdata[egdata$observer ==1, ], method = "ds", 
                      meta.data = list(width = 4))
   # run with Rsolnp, mono=TRUE, mono.strict=FALSE
   result.mono<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1), 
                    data = egdata[egdata$observer ==1, ], method = "ds", 
                    meta.data = list(width = 4,mono=TRUE,mono.strict=FALSE))
   # run with Rsolnp, mono=TRUE, mono.strict=FALSE
   result.strict<-ddf(dsmodel = ~mcds(key = "hn", formula = ~1), 
                      data = egdata[egdata$observer ==1, ], method = "ds", 
                      meta.data = list(width = 4,mono=TRUE,mono.strict=TRUE))
#     summary(result,se=TRUE)


   # test that the parameter is the same
   expect_that(result.optimx$par, equals(result.mono$par)) # FAIL 
   expect_that(result.optimx$par, equals(result.strict$par)) # FAIL
   expect_that(result.mono$par, equals(result.strict$par)) 

   # check that the summary gives the right answers too... 
   expect_that(summary(result.mono),prints_text("Monotonicity constraints were enforced."))
   expect_that(summary(result.strict),prints_text("Strict monotonicity constraints were enforced."))
   # check that the optimx() summary doesn't print anything about monotonicity
   expect_that(summary(result.optimx),prints_text("[^(Monotonicity constraints were enforced.)]"))
   expect_that(summary(result.optimx),prints_text("[^(Strict monotonicity constraints were enforced.)]"))
})

test_that("monotonicity warnings are correct", {
   # data setup  
   data(book.tee.data)
   region<<-book.tee.data$book.tee.region
   egdata<<-book.tee.data$book.tee.dataframe
   samples<<-book.tee.data$book.tee.samples
   obs<<-book.tee.data$book.tee.obs

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
   expect_that(summary(result),prints_text("[^(Monotonicity constraints were enforced.)]"))
   expect_that(summary(result),prints_text("[^(Strict monotonicity constraints were enforced.)]"))


})
