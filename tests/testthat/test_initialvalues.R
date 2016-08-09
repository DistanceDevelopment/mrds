library(mrds)

lnl.tol<-1e-3
par.tol<-1e-6

context("setting initial values")

test_that("initial values can be set",{

  data(book.tee.data)
  egdata<-book.tee.data$book.tee.dataframe

  # test that the initialvalues can be set
  result.mcds <- ddf(dsmodel=~mcds(key="hn", formula=~size),
                     data=egdata[egdata$observer==1,], method="ds",
                     control=list(initial=list(scale=c(0.5984749, 0.0212563))),
                     meta.data=list(width=4))
  expect_that(result.mcds$par, equals(setNames(c(0.5984749, 0.0212563),
                                               c("X.Intercept.", "size")),
                                      tolerance=1e-3))

  # same for hazard
  result.mcds <- ddf(dsmodel=~mcds(key="hr", formula=~size),
                     data=egdata[egdata$observer==1,], method="ds",
                     control=list(initial=list(scale=c(0.5984749, 0.0212563))),
                     meta.data=list(width=4))
  expect_that(result.mcds$par, equals(setNames(c(0.91510708, 0.66083586,
                                                 0.01749002),
                                               c("V1", "X.Intercept.", "size")),
                                      tolerance=1e-3))

})

