lnl.tol<-1e-3#4
par.tol<-1e-6

test_that("errors thrown for invalid models",{
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

  # errors when more pars than data
  egdata <- create.bins(egdata, 0:4)
  expect_error(ddf(dsmodel = ~cds(key = "hr", formula = ~1, adj.series = "cos",
                   adj.order = c(2, 3, 4), adj.scale = "width"), data = egdata,
                   method = "ds", meta.data = list(binned=TRUE,
                   breaks=0:4)),
               "More parameters to estimate than distance bins")

  obs <- data.frame(distance=rep(0, 17), object=1:17, observer=1)
  expect_error(a<-ddf(dsmodel = ~mcds(key = "hr", formula = ~1),
                   data = obs, method = "ds", meta.data = list(width=2)),
               "More parameters to estimate than unique distances")

})

