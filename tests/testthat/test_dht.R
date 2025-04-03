lnl.tol<-1e-3
par.tol<-1e-6

context("dht")

data(book.tee.data)
egdata <- book.tee.data$book.tee.dataframe
tee.data<-book.tee.data$book.tee.dataframe[book.tee.data$book.tee.dataframe$observer==1,]
ds.model <- ddf(dsmodel=~cds(key="hn", formula=~1), data=tee.data, method="ds",
                meta.data=list(width=4))

# same model, but calculating abundance
# need to supply the region, sample and observation tables
region <- book.tee.data$book.tee.region
samples <- book.tee.data$book.tee.samples
obs <- book.tee.data$book.tee.obs

test_that("golf tees", {

  # check that errors are thrown when the wrong ER variance is asked for
  expect_error(dht(ds.model, region.table=region,
                   sample.table=samples, obs.table=obs,
                   options=list(ervar="P2")),
               "Encounter rate variance estimator P2 / P3 may only be used with point transects, set with options=list(ervar=...)", fixed=TRUE)
  
  expect_error(dht(ds.model, region.table=region,
                   sample.table=samples, obs.table=obs,
                   options=list(ervar="P3")),
               "Encounter rate variance estimator P2 / P3 may only be used with point transects, set with options=list(ervar=...)", fixed=TRUE)
})


test_that("ptdata", {
  # fake up some pt data
  pt.sample <- data.frame(Sample.Label=1, Region.Label=1, Effort=1)
  data(ptdata.single)
  pt.obs <- data.frame(object=ptdata.single$object, Region.Label=1, Sample.Label=1)
  pt.region <- data.frame(Region.Label=1, Area=1)
  ptdata.single$distbegin <- (as.numeric(cut(ptdata.single$distance,10*(0:10)))-1)*10
  ptdata.single$distend <- (as.numeric(cut(ptdata.single$distance,10*(0:10))))*10
  model <- ddf(data=ptdata.single, dsmodel=~cds(key="hn"),
               meta.data=list(point=TRUE,binned=TRUE,breaks=10*(0:10), width=100))

  expect_warning(dht(model, region.table=pt.region, sample.table=pt.sample,
                     obs.table=pt.obs, options=list(ervar="O1")),
                 "Point transect encounter rate variance can only use estimators P2 or P3, switching to P2.", fixed=TRUE)

})


test_that("areas.supplied", {
  samples$CoveredArea<- samples$Effort * 2 * 4
  expect_equal(dht(ds.model, region, samples, obs,
                   options=list(areas.supplied=TRUE)),
               dht(ds.model, region, samples, obs))
})


test_that("dht with various opts after unif fixes", {
  
  # fit an independent observer model with point independence
  result.io <- ddf(dsmodel=~cds(key = "hn"), mrmodel=~glm(~distance),
                   data=egdata, method="io", meta.data=list(width=4))
  
  dht.result.io <- dht(result.io, 
                    region.table = region,
                    sample.table = samples,
                    obs.table = obs,
                    se = TRUE)
  
  expect_true(inherits(dht.result.io, "dht"))
  
  # fit an independent observer model with full independence
  result.io.fi <- ddf(mrmodel=~glm(~distance), data=egdata, method="io.fi",
                      meta.data=list(width = 4))
  
  dht.result.io.fi <- dht(result.io, 
                    region.table = region,
                    sample.table = samples,
                    obs.table = obs,
                    se = TRUE)
  
  expect_true(inherits(dht.result.io.fi, "dht"))
  
})


test_that("warning when only one sample", {
  library(Distance)
  data(capercaillie)
  capercaillie$size <- NULL
  conversion.factor <- convert_units("meter", "kilometer", "hectare")
  
  caper.ddf <- ddf(dsmodel = ~cds(key = "hn", formula = ~1),
                   data = capercaillie,
                   meta.data = list(width = 60))
  
  dht.tables <- Distance::unflatten(capercaillie)
  
  caper.dht <- expect_warning(dht(caper.ddf, 
                                  region.table = dht.tables$region.table,
                                  sample.table = dht.tables$sample.table,
                                  obs.table = dht.tables$obs.table,
                                  options = list(convert.units = conversion.factor)), "Only one sample, assuming abundance in the covered region is Poisson. See help on dht.se for recommendations.")
  
  caper.dht <- expect_warning(dht(caper.ddf, 
                                  region.table = dht.tables$region.table,
                                  sample.table = dht.tables$sample.table,
                                  obs.table = dht.tables$obs.table,
                                  options = list(convert.units = conversion.factor, varflag = 1)), "Only one sample, assuming variance of n is Poisson. See help on dht.se for recommendations.")
  
  caper.dht <- expect_no_warning(dht(caper.ddf, 
                                     region.table = dht.tables$region.table,
                                     sample.table = dht.tables$sample.table,
                                     obs.table = dht.tables$obs.table,
                                     options = list(convert.units = conversion.factor, varflag = 0)))
  
  # Make up a dataset so that there are multiple strata and some have only one transect
  # New dataset
  caper <- capercaillie
  caper$Region.Label <- "Strata A"
  caper2.a <- data.frame(Sample.Label = 2, Effort = 120, distance = abs(rnorm(45,0,30)),
                         object = NA, detected = 1, observer = 1, Region.Label = "Strata B",
                         Area = 1472)
  caper2.b <- data.frame(Sample.Label = 3, Effort = 120, distance = abs(rnorm(45,0,30)),
                         object = NA, detected = 1, observer = 1, Region.Label = "Strata B",
                         Area = 1472)
  caper2.c <- data.frame(Sample.Label = 4, Effort = 120, distance = abs(rnorm(90,0,30)),
                         object = NA, detected = 1, observer = 1, Region.Label = "Strata C",
                         Area = 1472)
  caper <- rbind(caper, caper2.a, caper2.b, caper2.c)
  caper$object <- 1:nrow(caper)
  
  # Model multi-strata data
  caper2.ddf <- ddf(dsmodel = ~cds(key = "hn", formula = ~1),
                    data = caper,
                    meta.data = list(width = 60))
  dht.tables2 <- Distance::unflatten(caper)
  
  caper2.dht <- expect_warning(dht(caper.ddf, 
                                   region.table = dht.tables2$region.table,
                                   sample.table = dht.tables2$sample.table,
                                   obs.table = dht.tables2$obs.table,
                                   options = list(convert.units = conversion.factor)),
                               "Only one sample in the following strata: Strata A, Strata C. For these strata, it is assumed abundance in the covered region is Poisson. See help on dht.se.")
  
  caper2.dht <- expect_warning(dht(caper.ddf, 
                                   region.table = dht.tables2$region.table,
                                   sample.table = dht.tables2$sample.table,
                                   obs.table = dht.tables2$obs.table,
                                   options = list(convert.units = conversion.factor, varflag = 1)),
                               "Only one sample in the following strata: Strata A, Strata C. For these strata, it is assumed variance of n is Poisson. See help on dht.se.")
  
  caper2.dht <- expect_no_warning(dht(caper.ddf, 
                                      region.table = dht.tables2$region.table,
                                      sample.table = dht.tables2$sample.table,
                                      obs.table = dht.tables2$obs.table,
                                      options = list(convert.units = conversion.factor, varflag = 0)))
  
})
