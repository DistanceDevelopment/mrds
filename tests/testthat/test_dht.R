lnl.tol<-1e-3
par.tol<-1e-6

context("dht")

test_that("golf tees",{

  data(book.tee.data)
  tee.data<-book.tee.data$book.tee.dataframe[book.tee.data$book.tee.dataframe$observer==1,]
  ds.model <- ddf(dsmodel=~cds(key="hn", formula=~1), data=tee.data, method="ds",
                  meta.data=list(width=4))

  # same model, but calculating abundance
  # need to supply the region, sample and observation tables
  region <- book.tee.data$book.tee.region
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs


  # check that errors are thrown when the wrong ER variance is asked for
  expect_error(dht(ds.model, region.table=region,
                   sample.table=samples, obs.table=obs,
                   options=list(ervar="P3")),
               "Encounter rate variance estimator P3 may only be used with point transects, set with options=list(ervar=...)", fixed=TRUE)

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
                 "Point transect encounter rate variance can only use estimators P2 or P3, switching to P3.", fixed=TRUE)

})

