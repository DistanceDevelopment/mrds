# what are the variance contributions for each stratum?
# based on the dht2 version
dht_variance_contributions <- function(res, groups, innes){

  # extract the columns we want
  CV_cont <- data.frame(ER          = res$summary$cv.ER,
                        Detection   = res$average.p.se/
                                      res$average.p)

  # if we're not doing Innes then we need to include group size
  if(!innes & !is.null(groups$se.Expected.S)){
    CV_cont$Groups <- groups$se.Expected.S/groups$Expected.S
  }

  # get the total
  CV_cont$Total <- sqrt(rowSums(CV_cont^2))

  # make that into percentages
  CV_cont <- (CV_cont^2/CV_cont[["Total"]]^2)*100
  CV_cont[["Total"]] <- NULL

  # zero ER contributions if only one sample
  CV_cont$ER[res$k==1] <- 0

  # sort and name
  CV_cont <- cbind(Region=res$summary$Region, CV_cont)

  return(CV_cont)
}
