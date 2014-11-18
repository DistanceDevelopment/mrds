# testing for summary.* functions

library(testthat)
library(mrds)

tol <- 1e-6

context("Summary tests")

test_that("summary results are correct",{

  # run the non-expect_that() code then the commented section below
  #  to "update" results
  #pt.result.summ<-summary(pt.result)
  #tee.result.summ<-summary(tee.result)
  #tee.result.trial.summ<-summary(tee.result.trial)
  #tee.result.trial.fi.summ<-summary(tee.result.trial.fi)
  #tee.result.rem.summ<-summary(tee.result.rem)
  #save(pt.result.summ,tee.result.summ,tee.result.trial.summ,
  #     tee.result.trial.fi.summ,tee.result.rem.summ,
  #     file="mrds/inst/testData/summaryresults.rda")

  ### load the results (summary objects)
  res.filename<-system.file("testData/summaryresults.rda", package="mrds")
  load(res.filename)


  ### io method -- line transects
  # load some data as usual
  data(book.tee.data)
  egdata <- book.tee.data$book.tee.dataframe

  # fit the model
  tee.result <- ddf(dsmodel = ~cds(key = "hn"), mrmodel = ~glm(~distance),
                data = egdata, method = "io", meta.data = list(width = 4))

  expect_that(summary(tee.result)$mr.summary,equals(tee.result.summ$mr.summary))
  expect_that(summary(tee.result)$ds.summary,equals(tee.result.summ$ds.summary))
  expect_that(summary(tee.result),equals(tee.result.summ))

  ### trial.fi method -- line transect
  tee.result.trial.fi <- ddf(mrmodel=~glm(formula=~distance),
                             dsmodel = ~mcds(key = "hn", formula = ~sex),
                             data=egdata, method="trial.fi",
                             meta.data=list(width=4))
  expect_that(summary(tee.result.trial.fi),equals(tee.result.trial.fi.summ))

  #### trial method -- line transect
  tee.result.trial <- ddf(mrmodel=~glm(formula=~distance),
                          dsmodel = ~mcds(key = "hn", formula = ~sex),
                          data=egdata, method="trial",
                          meta.data=list(width=4))
  expect_that(summary(tee.result.trial),equals(tee.result.trial.summ))


  ### rem method -- line transect
  tee.result.rem <- ddf(mrmodel=~glm(formula=~distance),
                          dsmodel = ~mcds(key = "hn", formula = ~sex),
                          data=egdata, method="rem",
                          meta.data=list(width=4))
  expect_that(summary(tee.result.rem),equals(tee.result.rem.summ))



  ## io method -- point data
  data(ptdata.dual)
  ptdata.dual$distbegin <- (as.numeric(cut(ptdata.dual$distance,
                                           10*(0:10)))-1)*10
  ptdata.dual$distend <- (as.numeric(cut(ptdata.dual$distance,10*(0:10))))*10

  pt.result <- ddf(method="io", data=ptdata.dual, dsmodel=~cds(key="hn"),
                   mrmodel=~glm(formula=~distance*observer),
                   meta.data=list(point=TRUE, binned=TRUE,
                                  breaks=10*(0:10),width=100))


  expect_that(summary(pt.result)$mr.summary,equals(pt.result.summ$mr.summary,
              tol=tol))
  expect_that(summary(pt.result)$ds.summary,equals(pt.result.summ$ds.summary,
              tol=tol))
  expect_that(summary(pt.result),equals(pt.result.summ,
              tol=tol))


})
