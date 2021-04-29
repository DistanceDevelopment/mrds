# general bug tests go here, they should go in another file if
# possible though

# TODO
# * check that qqplot.dff() gives a list of 4 when adjustments are used
# * need to surpress plotting?


# missing duplicates
context("bugs")

test_that("missing duplicates are detected",{
  data(book.tee.data)
  egdata <- book.tee.data$book.tee.dataframe

  # lose one data point
  egdata1 <- egdata[-nrow(egdata),]
  expect_error(ddf(mrmodel=~glm(~distance), data=egdata1, method="io.fi",
                   meta.data=list(width = 4)),
               "number of records for primary observer not equal to number for secondary observer")

  # lose two data point
  egdata2 <- egdata[c(1,nrow(egdata)),]
  expect_error(ddf(mrmodel=~glm(~distance), data=egdata2, method="io.fi",
                   meta.data=list(width = 4)),
               "objects 1, 162 do not have records for both observers")


})
