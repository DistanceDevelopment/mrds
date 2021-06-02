# Monotonicity constraints test file
# test that using Rsolnp gives the same results
# as optimx() when it should!

# boring bookexamples-based testing...
data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe

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


