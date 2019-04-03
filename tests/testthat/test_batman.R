# NA NA NA NA NA NA NA NA NA
library(mrds)

context("covariate NAs")

data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
egdata <- egdata[egdata$observer==1, ]

# make a missing value
egdata$sex[12] <- NA

test_that("ds",{
  expect_error(ddf(dsmodel=~mcds(key="hn", formula=~sex),
                   data=egdata, method="ds",
                   meta.data=list(width=4)),
               "NA covariate values in the data, check your data.")
})

data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe

egdata$sex[12] <- NA

test_that("io",{
  expect_error(ddf(dsmodel=~cds(key = "hn"), mrmodel=~glm(~distance+sex),
                   data=egdata, method="io", meta.data=list(width=4)),
               "NA covariate values in the data, check your data.")
})

test_that("io.fi",{
  expect_error(ddf(mrmodel=~glm(~distance+sex), data=egdata, method="io.fi",
                   meta.data=list(width = 4)),
               "NA covariate values in the data, check your data.")
})
test_that("rem",{
  expect_error(ddf(dsmodel=~cds(key = "hn"), mrmodel=~glm(~distance+sex),
                   data=egdata, method="rem", meta.data=list(width=4)),
               "NA covariate values in the data, check your data.")
})

test_that("rem.fi",{
  expect_error(ddf(mrmodel=~glm(~distance+sex), data=egdata, method="rem.fi",
                   meta.data=list(width = 4)),
               "NA covariate values in the data, check your data.")
})

#test_that("trial",{
#  expect_error(ddf(dsmodel=~cds(key = "hn"), mrmodel=~glm(~distance+sex),
#                   data=egdata, method="trial", meta.data=list(width=4)),
#               "NA covariate values in the data, check your data.")
#})
#
#test_that("trial.fi",{
#  expect_error(ddf(mrmodel=~glm(~distance+sex), data=egdata, method="trial.fi",
#                   meta.data=list(width = 4)),
#               "NA covariate values in the data, check your data.")
#})
