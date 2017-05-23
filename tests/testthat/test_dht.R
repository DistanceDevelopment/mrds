library(mrds)

lnl.tol<-1e-3#4
par.tol<-1e-6

context("dht")

test_that("golf tees",{

  data(book.tee.data)
  tee.data<-book.tee.data$book.tee.dataframe[book.tee.data$book.tee.dataframe$observer==1,]
  ds.model <- ddf(dsmodel=~cds(key="hn", formula=~1), data=tee.data, method="ds",
                  meta.data=list(width=4))

  ## Not run:

  # same model, but calculating abundance
  # need to supply the region, sample and observation tables
  region <- book.tee.data$book.tee.region
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs

  ds.dht.model <- dht(ds.model, region.table=region,
                      sample.table=samples, obs.table=obs)

res<-"Summary statistics:
  Region Area CoveredArea Effort   n  k        ER      se.ER      cv.ER
1      1 1040        1040    130  72  6 0.5538462 0.02926903 0.05284685
2      2  640         640     80  52  5 0.6500000 0.08292740 0.12758061
3  Total 1680        1680    210 124 11 0.5904762 0.03884115 0.06577936

Abundance:
  Label  Estimate       se         cv       lcl      ucl        df
1     1 123.22977 11.75088 0.09535745 101.72724 149.2774 43.918782
2     2  88.99928 13.37273 0.15025666  62.88926 125.9495  7.658529
3 Total 212.22905 21.33325 0.10051992 173.30068 259.9019 40.063061

Density:
  Label  Estimate         se         cv        lcl       ucl        df
1     1 0.1184902 0.01129892 0.09535745 0.09781465 0.1435359 43.918782
2     2 0.1390614 0.02089490 0.15025666 0.09826447 0.1967961  7.658529
3 Total 0.1263268 0.01269836 0.10051992 0.10315516 0.1547035 40.063061

Summary for individuals

Summary statistics:
  Region Area CoveredArea Effort   n       ER     se.ER      cv.ER mean.size
1      1 1040        1040    130 229 1.761538 0.1165805 0.06618107  3.180556
2      2  640         640     80 152 1.900000 0.3342319 0.17591151  2.923077
3  Total 1680        1680    210 381 1.814286 0.1391400 0.07669132  3.072581
    se.mean
1 0.2086982
2 0.2261991
3 0.1537082

Abundance:
  Label Estimate       se        cv      lcl      ucl        df
1     1 391.9391 40.50494 0.1033450 317.2772 484.1706 27.423281
2     2 260.1517 50.20666 0.1929899 162.2494 417.1289  5.786774
3 Total 652.0909 73.79806 0.1131714 516.5938 823.1274 23.815561

Density:
  Label  Estimate         se        cv       lcl       ucl        df
1     1 0.3768645 0.03894706 0.1033450 0.3050742 0.4655487 27.423281
2     2 0.4064871 0.07844791 0.1929899 0.2535147 0.6517639  5.786774
3 Total 0.3881493 0.04392741 0.1131714 0.3074963 0.4899568 23.815561

Expected cluster size
  Region Expected.S se.Expected.S cv.Expected.S
1      1   3.180556     0.2114629    0.06648615
2      2   2.923077     0.1750319    0.05987935
3  Total   3.072581     0.1391365    0.04528327"

  expect_output(print(ds.dht.model), res)

})

