# testing for plotting functions

library(testthat)
library(mrds)

context("Plotting tests: average.line")

test_that("average.line for plot.io is correct (unconditional df)",{

  # load some data as usual
  data(book.tee.data)
  region <- book.tee.data$book.tee.region
  egdata <- book.tee.data$book.tee.dataframe
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs

  # fit the model
  result <- ddf(dsmodel = ~cds(key = "hn"), mrmodel = ~glm(~distance),
                data = egdata, method = "io", meta.data = list(width = 4))

  # make the grid (as in plot.io)
  divisions <- 25 # default in plot.io
  finebr <- (4/divisions)*(0:divisions)

  # make the line
  aline <- mrds:::average.line(finebr,obs=1,result)

  # "true" values
  aline.truth <- c(0.942984310350384,0.931758960736299,0.916033134386296,
                   0.895960681152277,0.871766127423631,0.843743459341331,
                   0.812252811904826,0.777714776947772,0.740602126839018,
                   0.701428901443524,0.660737024685639,0.619080889235704,
                   0.577010639634571,0.535055145122935,0.493705825337338,
                   0.453402523523292,0.414522484905807,0.377373199698558,
                   0.342189456128584,0.309134492201619,0.278304718215714,
                   0.249737174975323,0.223418735462668,0.199296054654547,
                   0.177285396274316)

  expect_that(aline$value,equals(aline.truth,tol=1e-8))

})
