lnl.tol<-1e-3#4
par.tol<-1e-6

context("errors thrown for invalid models")
test_that("golf tee data gives the same results as Distance",{
  data(book.tee.data)
  egdata<-book.tee.data$book.tee.dataframe


  # half-normal with order 1 cosine adjustment
  expect_error(result.cds<-ddf(dsmodel = ~cds(key = "hn",
                                              adj.series="cos",adj.order=1),
                               data = egdata, method = "ds",
                               meta.data = list(width = 4)),
               "Cosine adjustments must be of order >2 for half-normal key functions")

  # half-normal with order 1 hermite adjustment
  expect_error(result.cds<-ddf(dsmodel = ~cds(key = "hn",
                                              adj.series="herm",adj.order=1),
                               data = egdata, method = "ds",
                               meta.data = list(width = 4)),
               "Hermite polynomial adjustment terms of order < 4 selected")

  # hazard-rate with order 1 cosine adjustment
  expect_error(result.cds<-ddf(dsmodel = ~cds(key = "hr",
                                              adj.series="cos",adj.order=1),
                               data = egdata, method = "ds",
                               meta.data = list(width = 4)),
               "Cosine adjustments must be of order >2 for hazard-rate key functions")

  # hazard-rate with order 1 simple poly adjustment
  expect_error(result.cds<-ddf(dsmodel = ~cds(key = "hr",
                                              adj.series="poly",adj.order=1),
                               data = egdata, method = "ds",
                               meta.data = list(width = 4)),
               "Odd polynomial adjustment terms selected")

})

