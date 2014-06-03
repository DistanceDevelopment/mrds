# bugs in trial and trial.io go here...

library(testthat)
library(mrds)

context("trial and trial.fi")

# bug from Louise Burt 31/10/11
test_that("detected and missed are the same in summary and det.tables", {

   # bug appeared because there were zeros in the data but when using cut()
   # you must specify include.lowest=TRUE for them to be included.

   # grab some data
   data(book.tee.data)
   egdata<<-book.tee.data$book.tee.dataframe


   # no zeros, just check it works anyway
   xx<-ddf(mrmodel=~glm(formula=~distance),
           dsmodel = ~mcds(key = "hn", formula = ~sex),
           data = egdata, method = "trial.fi", meta.data = list(width = 4))

   tabs<-det.tables(xx,breaks=c(0,.5,1,2,3,4))
   summ<-summary(xx)

   # the sums of the columns in tabs should match the summary numbers
   expect_that(sum(tabs$Observer1[,2]), equals(summ$n1))
   expect_that(sum(tabs$Observer2[,2]), equals(summ$n2))
   expect_that(sum(tabs$Duplicates), equals(summ$n3))


   # there are no zeros, so fix that and re-run the model
   egdata$distance[egdata$distance==0.02]<-0
   xx<-ddf(mrmodel=~glm(formula=~distance),
           dsmodel = ~mcds(key = "hn", formula = ~sex),
           data = egdata, method = "trial.fi", meta.data = list(width = 4))
   tabs<-det.tables(xx,breaks=c(0,.5,1,2,3,4))
   summ<-summary(xx)
   # test the bug
   expect_that(sum(tabs$Observer1[,2]), equals(summ$n1))
   expect_that(sum(tabs$Observer2[,2]), equals(summ$n2))
   expect_that(sum(tabs$Duplicates), equals(summ$n3))


})

# test results from golf tee data
test_that("results from golf tee data (trial.fi mode) work (distance covar)", {

  data(book.tee.data)
  # extract the list elements from the dataset into easy-to-use objects
  detections <- book.tee.data$book.tee.dataframe
  # make sure sex and exposure are factor variables
  detections$sex <- as.factor(detections$sex)
  detections$exposure <- as.factor(detections$exposure)
  region <- book.tee.data$book.tee.region
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs

  # Fit the model
  fi.mr.dist <- ddf(method="trial.fi",
                    mrmodel=~glm(link="logit",formula =~distance),
                    data = detections, meta.data = list(width = 4))

  # setup the "true" parameter estimates
  par.ests <- c(2.900233,-1.058677)
  names(par.ests) <- c("(Intercept)","distance")

  # test parameters
  expect_that(fi.mr.dist$par,equals(par.ests,tolerance=1e-6))

  # test likelihood
  expect_that(fi.mr.dist$lnl,equals(-224.4047,tolerance=1e-6))

  # test average p
  expect_that(summary(fi.mr.dist)$average.p,equals(0.6423252,tolerance=1e-6))

})

test_that("results from golf tee data (trial.fi mode) work (distance+sex*exposure covar)", {
  data(book.tee.data)
  # extract the list elements from the dataset into easy-to-use objects
  detections <- book.tee.data$book.tee.dataframe
  # make sure sex and exposure are factor variables
  detections$sex <- as.factor(detections$sex)
  detections$exposure <- as.factor(detections$exposure)
  region <- book.tee.data$book.tee.region
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs

  fi.mr.dist.sex.exposure <- ddf(method = "trial.fi",
                               mrmodel=~glm(link = "logit",
                               formula = ~distance + sex * exposure),
                               data = detections, meta.data = list(width = 4))

  par.ests <- c(-0.1938752,-1.9137954,3.9599369,4.7888851,-1.7273873)
  names(par.ests) <- c("(Intercept)","distance","sex1","exposure1",
                       "sex1:exposure1")

  # test parameters
  expect_that(fi.mr.dist.sex.exposure$par,equals(par.ests,tolerance=1e-6))

  # test likelihood
  expect_that(fi.mr.dist.sex.exposure$lnl,equals(-196.9022,tolerance=1e-6))

  # test average p
  expect_that(summary(fi.mr.dist.sex.exposure)$average.p,
              equals(0.5245967,tolerance=1e-6))


})
