# bugs in trial and trial.io go here...

library(testthat)
library(mrds)

context("trial and trial.io")

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
