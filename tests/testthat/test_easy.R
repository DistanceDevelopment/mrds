library(mrds)

lnl.tol<-1e-3#4
par.tol<-1e-6

context("easy tests: golf tees")
# taken from Laura Marshall's RUnit tests

test_that("golf tee data gives the same results as Distance",{

  data(book.tee.data)
  egdata<-book.tee.data$book.tee.dataframe


  #CDS
  result.cds <- ddf(dsmodel=~cds(key="hn"), mrmodel=~glm(~distance),
                    data=egdata, method="io", meta.data=list(width=4))
  expect_that(result.cds$Nhat, equals(232.0015, tolerance=1e-6))

  #GLM
  result.glm <- ddf(mrmodel=~glm(~distance), data=egdata, method="io.fi",
                  meta.data=list(width=4))
  expect_that(result.glm$Nhat, equals(186.0947, tolerance=1e-6))

  #MCDS
  #checkException(ddf(dsmodel = ~mcds(key = "hn", formula = ~1), data = newdata, method = "ds", meta.data = list(width = 4)))

  result.mcds <- ddf(dsmodel=~mcds(key="hn", formula=~1),
                   data=egdata[egdata$observer==1, ], method="ds",
                   meta.data=list(width=4))
  expect_that(result.mcds$Nhat, equals(212.229, tolerance=1e-3))

  #rr <- ddf(dsmodel = ~mcds(key = "hn", formula = ~sex),
  #          data=subset(egdata, detected==1 & observer==1), method="ds",
  #          meta.data=list(width=4))

#  # check that uniform key works
#  result.unif<-ddf(dsmodel = ~cds(key = "unif",adj.series="cos",adj.order=2), 
#                  data = egdata, method = "ds", meta.data = list(width = 4))
#  expect_that(result.unif$par, equals(0.7005934,tolerance=1e-6))
  
# there should be an error if we don't supply adjustments with uniform key
# This was removed because it now works
#  expect_error(ddf(dsmodel = ~cds(key = "unif"), 
#                  data = egdata, method = "ds", meta.data = list(width = 4)))

})

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

context("easy tests: line transect example")
# line transect example from Distance
### check likelihood and pars
# some of these give better likelihoods (lower -logL) than
# Distance, so just test that...

# code here is a bit cryptic...
# only want this data visible from the tests, so need to load it this way...
# load the data
ex.filename<-system.file("testData/ltexample.rda", package="mrds")
load(ex.filename)
# and the model definitions and results from distance
res.filename<-system.file("testData/ltresults.rda", package="mrds")
load(res.filename)

width<-25

# load the results from Distance
# same order as here, so just iterate over the list

test.df<-function(mcds.bit,dat,width,mono=FALSE,strict=FALSE,mono.tol=1e-7,showit=0){
  suppressMessages(ddf(dsmodel=mcds.bit,data=dat,method="ds",
      meta.data=list(width=width,mono=mono,mono.strict=strict),
#      control=list(mono.tol=mono.tol,showit=showit))
      control=list(mono.tol=1e-5,mono.delta=1e-5,showit=showit)))
}

models<-ltmodels
model.set<-1:nrow(models)

### weird bugs
#wb<-c(11,12,27,28,30)
### really weird bugs
#rwb<-c(21)

# remove some models
# 12 has a paramter ~=0, so don't try to fit that
# 30 is really weird:
#  ack<-function(x){exp(-x^2/(2*2.5^2))*(1+-0.9*((x/2.5)^4))}
#  plot(seq(0,25,len=1000),ack(seq(0,25,len=1000)))
#  integrate(ack,upper=25,lower=0)

model.set<-model.set[-c(1,10,11,12,14,16,21,22,27:30)]

better<-rep(0,nrow(models))

for(i in model.set){

  set.seed(1245)
  # uncomment for debug
  #cat("\nmodel",i,":\n")

  # construct the model
  mcds.call<-paste("~mcds(key=\"")
  key<-switch(models$key[i],HN="hn",HA="hr")
  mcds.call<-paste(mcds.call,key,"\"",sep="")
  adjust<-switch(models$adj[i],
                 CO="cos",
                 PO="poly",
                 HE="herm")
  scaling<-models$scaling[i]
  if(scaling=="sigma") scaling<-"scale"
  # build the bits relating to adjustments
  if(models$noAdj[i]!=0){
    adj.order<-models$order[i]
    if(grepl(",",adj.order)){
      adj.order<-paste("c(",adj.order,")",sep="")
    }
    mcds.call<-paste(mcds.call,",adj.series=\"",adjust,"\",
                                 adj.order=",adj.order,",
                                 adj.scale=\"",scaling,"\"",sep="")
  }
  # monotonicity bits
  mono<-switch(models$monotone[i],
               None=FALSE,
               Weak=TRUE,
               strict=TRUE)
  mono.strict<-switch(models$monotone[i],
                None=FALSE,
                Weak=FALSE,
                strict=TRUE)

  mcds.call<-paste(mcds.call,",formula=~1)",sep="")

  this.model<-ltresults[[i]]

  # actually fit some models
  if(this.model$status==1 | this.model$status==2){
    result<-try(test.df(eval(parse(text=mcds.call)),ltexample,width,
                    mono=mono,strict=mono.strict,showit=0))

    expect_that(all(class(result)=="try-error"),is_false(), label=i)

    if(all(class(result)!="try-error")){

      this.test<-paste(mcds.call,"\nmono=",mono,
                                 "\nmono.strict=",mono.strict,"\n")

      if(result$lnl<=this.model$LnL){
        test_that(this.test,{
          expect_that(result$lnl<=this.model$LnL,is_true(),
                      label=paste("Likelihood for model",i,
                                  "better than MCDS"))
        })
      }else{
        test_that(this.test,{
          expect_that(result$lnl,equals(this.model$LnL,tol=lnl.tol),
                      label=paste("Likelihood for model",i,
                                  "the same as MCDS"))
        })
      }

      #if(result$lnl<=this.model$LnL){
      # # make a note of when mrds does better than Distance
      # better[i]<-1

      #}else{
      #  test_that(paste(mcds.call,"pars"),{
      #    names(result$par)<-NULL
      #    # do something about the first one/two pars of hn/hr
      #    this.model$param[1]<-log(this.model$param[1])
      #    if(key=="hr") this.model$param[2]<-log(this.model$param[2])

      #    expect_that(result$par, equals(this.model$param,tolerance=par.tol))
      #  })
      #}
    }

  # uncomment for debug
  #cat("\n")

  }else{
    cat(paste("\n",i," -- MCDS error!\n"))
  }
}

# extra column of when mrds beat Distance
#models<-cbind(models,better)

#cat("the following were fine:\n")
#print(models[model.set,])
#
#cat("weird bugs:\n")
#print(models[wb,])
#
#cat("not as good likelihood:\n")
#print(models[fails,])

#
#
##context("easy tests: ducknest")
#  dat<-c(0.06,0.07,0.04,0.01,0.37,0.36,0.51,0.45,0.32,0.61,0.61,
#         0.66,0.69,1.02,1.15,1,1.03,1.05,1.41,1.4,1.69,1.61,2,1.97,
#         1.95,2.13,2.27,0.28,0.23,0.1,0.38,0.59,0.86,0.86,0.89,0.82,
#         0.76,1.19,1.05,0.94,0.94,1.27,1.23,1.41,1.39,1.31,1.74,2.21,
#         2.29,2.19,2.38,2.22,0.2,0.1,0.06,0.21,0.13,0.23,0.33,0.39,
#         0.35,0.44,0.42,0.38,0.74,0.66,0.68,0.8,0.7,0.76,0.89,1.07,
#         0.96,0.98,0.92,1.26,1.47,1.31,1.22,1.67,1.81,2.33,2.23,2.14,
#         0.17,0.19,0.25,0.03,0.13,0.18,0.57,0.53,0.51,0.52,0.42,0.55,
#         0.63,0.89,0.75,0.77,0.76,0.88,1.09,1.09,1.21,1.4,1.58,1.61,
#         1.64,1.66,1.65,1.95,1.91,2.3,2.11,2.18,0.06,0.13,0.12,0.49,
#         0.48,0.54,0.46,0.35,0.53,0.74,1.07,1.1,1.02,1.34,1.4,1.29,
#         1.43,1.37,1.38,1.23,1.69,1.68,1.7,1.61,1.56,1.86,1.83,0.1,
#         0.08,0.03,0.18,0.24,0.42,0.52,0.41,0.34,0.63,0.71,0.75,0.75,
#         0.97,1.08,1.12,1.31,1.33,1.4,1.22,1.71,1.59,1.58,2.05,1.91,
#         1.84,1.99,2.32,2.36,2.24,2.12,0.07,0.27,0.09,0.05,0.03,0.24,
#         0.44,0.32,0.42,0.33,0.73,0.78,0.79,0.62,0.74,0.87,0.67,1.17,
#         1.08,1.49,1.42,1.34,1.68,1.58,1.67,1.94,2,2.09,2.37,2.16,
#         0.08,0.16,0.21,0.09,0.51,0.56,0.75,0.62,1.1,1.17,1.16,1.12,
#         1.43,1.42,1.71,1.71,1.64,1.66,1.82,2.08,1.88,2.17,2.35,2.39,
#         0.05,0.03,0.51,0.38,0.58,0.45,0.32,0.72,1.08,1.13,1.12,1.18,
#         1.01,1.05,1.38,1.39,1.25,1.31,1.68,1.73,1.81,1.98,1.95,2.09,
#         2.15,2.27,0.19,0.23,0.25,0.31,0.68,1.02,0.97,1.15,0.91,1.41,
#         1.21,1.42,1.7,1.66,1.61,1.89,2.16,0.17,0.48,0.59,0.58,0.59,
#         0.63,0.85,0.64,0.67,1,1.41,1.38,1.26,1.35,1.78,1.55,1.58,
#         1.73,1.72,1.86,2.06,2,2.32,2.2,2.25,0.27,0.31,0.77,0.83,
#         0.72,0.64,1.02,1.51,1.76,1.69,2,1.91,1.94,2.07,1.93,1.9,
#         1.88,2.27,2.25,2.36,0.21,0.26,0.05,0.53,0.34,0.82,0.61,0.79,
#         0.86,0.77,0.69,0.74,0.85,1.71,1.68,1.64,1.85,1.98,1.88,2.11,
#         2.18,2.4,0.07,0.04,0.45,0.49,0.54,0.57,0.41,0.53,0.54,0.4,
#         0.66,0.86,1.16,1.18,1.05,0.96,0.91,1.35,1.4,1.46,1.3,1.29,
#         1.68,1.81,2.35,2.29,0.27,0.09,0.2,0.21,0.35,0.57,0.48,0.47,
#         0.69,0.7,0.67,0.68,0.97,0.96,0.98,1.14,0.92,1.02,1.37,1.43,
#         1.4,1.22,1.26,1.44,1.54,1.67,2.07,1.86,1.95,2,2.05,2.12,
#         0.28,0.12,1.02,1.1,1.04,1.22,1.39,1.28,1.23,1.26,1.39,1.48,
#         1.39,1.7,1.9,1.99,1.87,2.31,2.27,0.22,0.18,0.22,0.26,0.36,
#         0.59,0.34,0.36,0.86,0.69,0.67,0.98,0.97,1.06,1.19,1.14,1.45,
#         1.35,1.32,1.48,1.78,1.69,1.62,1.55,1.66,1.71,2.05,2.06,2.38,
#         2.17,2.38,0.24,0.06,0.23,0.23,0.13,0.28,0.37,0.39,0.32,0.47,
#         0.36,0.44,0.84,0.74,0.66,0.85,1.31,1.43,1.33,1.39,1.3,1.36,
#         1.75,1.72,2,2.19,2.38,2.27,2.37,0.23,0.26,0.28,0.03,0.33,
#         0.49,0.7,0.61,0.76,0.7,0.86,1.14,1.05,1.01,1.31,1.49,1.37,
#         1.29,1.31,1.52,1.61,1.51,1.76,1.63,1.84,2.29,2.27,0.19,0.23,
#         0.27,0.08,0.19,0.55,0.62,0.71,0.73,0.76,0.85,0.84,0.71,1.06,
#         1.18,1,0.94,1.05,1.31,1.38,1.27,1.21,1.51,1.84,1.92,1.94,
#         2.18,2.26,2.36,2.32,2.38,2.13)
#
#  dat<-data.frame(distance=dat,object=1:length(dat),
#                  observed=rep(1,length(dat)))
#width<-max(dat$distance)
#
##
##test_that("ducknest gives the right answers",{
##
#test_that("half-normal",{
#
#  result<-test.df(~mcds(key="hn",formula=~1),dat,width)
#  expect_that(result$lnl<=-317.49061,is_true())
#  #expect_that(result$par, equals(2.300547,tolerance=1e-6))
#
#  result<-test.df(~mcds(key="hn",formula=~1),dat,width,mono=TRUE,strict=FALSE)
#  expect_that(result$lnl<=-317.49061,is_true())
#  #expect_that(result$par, equals(2.300547,tolerance=1e-6))
#
#  result<-test.df(~mcds(key="hn",formula=~1),dat,width,mono=TRUE,strict=TRUE)
#  expect_that(result$lnl<=-317.49061,is_true())
#  #expect_that(result$par, equals(2.300547,tolerance=1e-6))
#})
#
##})
