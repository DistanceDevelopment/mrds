# hopefully these should cause some problems
library(testthat)
library(mrds)
# pull in the mixture simulation code

context("nasty mixtures")
test_that("easy mixture",{

   # save the seed
   #rngsave<-.Random.seed
   set.seed(3141)

   # simulate some data
##!##   dat<-sim.mix(500,c(0.7,0.05),c(0.6,0.4),4,means=c(0,0))
##!##   dat<-data.frame(distance=dat,object=1:length(dat),observed=rep(1,length(dat)))
##!##
##!##   trunc<-2
##!##
##!##   # fit without constraint
##!##   result.n<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=FALSE))
##!##
##!##   result.s<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=TRUE))
##!##
##!##   result.w<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
})

test_that("less easy mixture -- mode away from zero",{

   # save the seed
   #rngsave<-.Random.seed
   set.seed(3141)
   # simulate some data
##!##   dat<-sim.mix(500,c(0.3,1),c(0.6,0.4),10,means=c(3,0))
##!##
##!##
##!##   dat<-data.frame(distance=dat,object=1:length(dat),observed=rep(1,length(dat)))
##!##
##!##   trunc<-4.5
##!##
##!##   # fit cosine adjsutment 
##!##   result.n<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=FALSE))
##!##
##!##   result.s<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=TRUE))
##!##
##!##   result.w<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
##!##
##!##
##!##
##!##   # hazard rate
##!##   result.n<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=FALSE))
##!##
##!##   result.s<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=TRUE))
##!##
##!##   result.w<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
##!##                 data=dat, method = "ds", 
##!##                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
})

#test_that("less easy mixture -- spike",{
#
#   # save the seed
#   #rngsave<-.Random.seed
#   set.seed(3141)
#   # simulate some data
#   dat<-sim.mix(500,c(0.3,2),c(0.9,0.1),10,means=c(0,0))
#
#   dat<-data.frame(distance=dat,object=1:length(dat),observed=rep(1,length(dat)))
#
#   trunc<-4.5
#
#   # fit cosine adjsutment 
#   result.n<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=FALSE))
#
#   result.s<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=TRUE))
#
#   result.w<-ddf(dsmodel = ~mcds(key = "hn",formula=~1,adj.series="cos",adj.order=c(2),adj.scale="scale"), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
#
#
## hazard rate
#   result.n<-ddf(dsmodel=~mcds(key="hr",formula=~1),data=dat, method="ds", 
#                 meta.data=list(width=trunc,mono=FALSE))
#
#   result.w<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
#
#   result.s<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=TRUE,mono.strict=TRUE))
#par(mfrow=c(1,3))
#plot(result.n)
#plot(result.s)
#plot(result.w)
#
#   # make it more difficult -- this all breaks
#   dat<-c(dat$distance,rep(0,140))
#   dat<-data.frame(distance=dat,object=1:length(dat),observed=rep(1,length(dat)))
#
#
## hazard rate
#   result.n<-ddf(dsmodel=~mcds(key="hr",formula=~1),data=dat, method="ds", 
#                 meta.data=list(width=trunc,mono=FALSE))
#
#   result.w<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))
#
#   result.s<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
#                 data=dat, method = "ds", 
#                 meta.data=list(width=trunc,mono=TRUE,mono.strict=TRUE))
#par(mfrow=c(1,3))
#plot(result.n)
#plot(result.s)
#plot(result.w)
#
#
#})

#   # do something weird -- this messes everything up

   # this breaks
   ## Jeff's suggestion, generate data then add a whole bunch of zeros... 
   #dat<-sim.mix(40,c(1),c(0.6),4,means=c(0,0))
   #dat<-c(dat,rep(0,15))
   #dat<-data.frame(distance=dat,object=1:length(dat),observed=rep(1,length(dat)))

   #trunc<-2.5
   #result.n<-ddf(dsmodel=~mcds(key="hr",formula=~1),data=dat, method="ds", 
   #              meta.data=list(width=trunc,mono=FALSE),control=list(lowerbounds=c(0.0001,-2.5)))

   #result.n<-ddf(dsmodel=~mcds(key="hr",formula=~1),data=dat, method="ds", 
   #              meta.data=list(width=trunc,mono=FALSE),
   #              control=list(lowerbounds=c(result.n$par[1]-1e-4,-10),upperbounds=c(result.n$par[1]+1e-4,10),initial=list(scale=result.n$ds$aux$ddfobj$scale$parameters, shape=result.n$ds$aux$ddfobj$shape$parameters)))


   #result.w<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
   #              data=dat, method = "ds", 
   #              meta.data=list(width=trunc,mono=TRUE,mono.strict=FALSE))

   #result.s<-ddf(dsmodel = ~mcds(key = "hr",formula=~1), 
   #              data=dat, method = "ds", 
   #              meta.data=list(width=trunc,mono=TRUE,mono.strict=TRUE))

   # recover the seed
   #.Random.seed<-rngsave



