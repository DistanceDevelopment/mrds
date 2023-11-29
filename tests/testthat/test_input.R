# Input checks
context("user input")

test_that("MCDS optimizer is not used with double observer analyses",{
  
  data(book.tee.data)
  region <- book.tee.data$book.tee.region
  egdata <- book.tee.data$book.tee.dataframe
  samples <- book.tee.data$book.tee.samples
  obs <- book.tee.data$book.tee.obs
  
  expect_warning(ddf(dsmodel=~cds(key = "hn"), mrmodel=~glm(~distance),
                   data=egdata, method="io", meta.data=list(width=4), 
                   control = list(optimizer = "MCDS")),
                 "The MCDS optimizer cannot currently be used with double observer analyses, the R opitimizer will be used instead.")
  
})
