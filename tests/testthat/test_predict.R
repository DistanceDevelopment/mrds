# test predictions with covariates


data(book.tee.data)

context("prediction with covars")


# fit without covar
dd <- ddf(~mcds(formula=~1, key="hn"), meta.data=list(width=4),
          data=book.tee.data$book.tee.dataframe)

# fit with covar
dd_h <- ddf(~mcds(formula=~as.factor(sex), key="hn"), meta.data=list(width=4),
          data=book.tee.data$book.tee.dataframe)

nd <- data.frame(distance=0)
test_that("simple prediction, no covars",{

  pp <- predict(dd, newdata=nd)

  expect_equal(pp$fitted, unname(fitted(dd)[1]))
})

test_that("covar df, no covar in data",{
  expect_error(predict(dd_h, newdata=nd), "columns in `newdata` do not match those in fitted model")
})



nd <- data.frame(distance=0,
                 sex = book.tee.data$book.tee.dataframe[1,]$sex)

test_that("covar df, covar in data",{
  pp <- predict(dd_h, newdata=nd)
  expect_equal(pp$fitted, unname(fitted(dd_h)[1]))
})



# test that a factor outside of the observed ones makes things explode

nd <- data.frame(distance=0,
                 sex = 6)
test_that("covar df, no covar level in data",{
  expect_error(predict(dd_h, newdata=nd), "fields or factor levels in `newdata` do not match data used in fitted model")
})

