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

# make prediction where there are NAs in sex
nd <- data.frame(distance=c(0,0),
                 sex = c(book.tee.data$book.tee.dataframe[1,]$sex, NA))

test_that("covar df, NA covar in data",{
  pp <- predict(dd_h, newdata=nd)
  expect_equal(pp$fitted, c(unname(fitted(dd_h)[1]), NA))
})

# make prediction where there are NAs in sex
nd <- data.frame(distance=0,
                 sex = c(book.tee.data$book.tee.dataframe$sex[1], NA))

test_that("covar df, NA covar in data",{
  pp <- predict(dd_h, newdata=nd)
  expect_equal(pp$fitted, c(unname(fitted(dd_h)[1]), NA))
})

# extra nonsense, need to ignore additional columns not in df
nd <- data.frame(distance=c(0,0),
                 sex = book.tee.data$book.tee.dataframe[1:2,]$sex,
                 foo = c(NA, NA))
test_that("covar df, NA in noncovar data",{
  pp <- predict(dd_h, newdata=nd)
  expect_equal(pp$fitted, c(unname(fitted(dd_h)[1:2])))
})

test_that("create.bins cuts data correctly - called from predict",{
  # Testing fix for issue #73
  data(book.tee.data)
  tee.data <- subset(book.tee.data$book.tee.dataframe, observer==1)
  tee.data$distance <- 1000*tee.data$distance
  cutpoints = c(0,250,500,750,1000, 1500, 2000, 3000, 4000)
  dat.mrds <- create.bins(tee.data, cutpoints = cutpoints)
  expect_false(any(is.na(dat.mrds$distbegin)))
  expect_false(any(is.na(dat.mrds$distend)))
})

test_that("Predict for hr model with covariates when se.fit = TRUE",{
  # Testing fix for issue #84
  dat <- data.frame(distance=abs(rnorm(100)), x=rnorm(100))
  hrcov <- Distance::ds(dat, formula=~x, key="hr")
  nd <- data.frame(distance = 0.25, x=0.5)
  tmp <- predict(hrcov, nd, se.fit=TRUE)
  expect_equal(length(tmp), 2)
})

