#' Quantile-quantile plot and goodness of fit tests for detection functions
#'
#' Constructs a quantile-quantile (Q-Q) plot for fitted model as a graphical check of goodness of fit. Formal goodness of fit testing for detection function models using Kolmogorov-Smirnov and Cramer-von Mises tests. Both tests are based on looking at the quantile-quantile plot produced by \code{\link{qqplot.ddf}} and deviations from the line x=y.
#'
#' The Kolmogorov-Smirnov test asks the question "what's the largest vertical distance between a point and the y=x line?" It uses this distance as a statistic to test the null hypothesis that the samples (EDF and CDF in our case) are from the same distribution (and hence our model fits well). If the deviation between the y=x line and the points is too large we reject the null hypothesis and say the model doesn't have a good fit.
#'
#' Rather than looking at the single biggest difference between the y=x line and the points in the Q-Q plot, we might prefer to think about all the differences between line and points, since there may be many smaller differences that we want to take into account rather than looking for one large deviation. Its null hypothesis is the same, but the statistic it uses is the sum of the deviations from each of the point to the line.
#
#' @section Details:

#' Note that a bootstrap procedure is required to ensure that the p-values from the procedure are correct as the we are comparing the cumulative distribution function (CDF) and empirical distribution function (EDF) and we have estimated the parameters of the detection function.

#'
#' @param model fitted distance detection function model object
#' @param plot the Q-Q plot be plotted or just report statistics?
#' @param nboot number of replicates to use to calculate p-values for the goodness of fit test statistics
#' @param ks perform the Kolmogorov-Smirnov test (this involves many bootstraps so can take a while)
#' @param \dots additional arguments passed to \code{\link{plot}}
#' @export
#' @return A list of goodness of fit related values: \item{edf}{matrix of lower
#'   and upper empirical distribution function values} \item{cdf}{fitted
#'   cumulative distribution function values} \item{ks}{list with K-S statistic
#'   (\code{Dn}) and p-value (\code{p})} \item{CvM}{list with CvM statistic
#'   (\code{W}) and p-value (\code{p})}
#' @author Jeff Laake, David L Miller
#' @seealso \code{\link{ddf.gof}}, \code{\link{cdf.ds}}
#' @references Burnham, K.P., S.T. Buckland, J.L. Laake, D.L. Borchers, T.A.
#'   Marques, J.R.B. Bishop, and L. Thomas. 2004.  Further topics in distance
#'   sampling. pp: 385-389. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'   D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'   Oxford University Press.
#' @keywords utility
#' @importFrom graphics abline
qqplot.ddf <- function(model, plot=TRUE, nboot=100, ks=FALSE, ...){

  # get the edf/cdf values
  edf_cdf <- get_edf_cdf(model)

  if(plot){
    plot(edf_cdf$upper.edf, edf_cdf$cdfvalues,
         xlab="Empirical cdf", ylab="Fitted cdf",
         xlim=c(0,1), ylim=c(0,1), asp=1, ...)
    abline(0, 1, ...)
  }

  # don't do KS if the model is not "ds"
  if(!("ds" %in% class(model)) & ks){
    warning("Can't calculate Kolmogorov-Smirnov p-value for non-\"ds\" models, only calculating Cramer-von Mises")
    ks <- FALSE
  }

  # do the tests
  gof_p <- gof.pvalues(model, ks=ks, nboot=nboot)

  # build return object
  return(list(edf=cbind(edf_cdf$lower.edf, edf_cdf$upper.edf),
              cdf=edf_cdf$cdfvalues,
              ks=list(Dn=gof_p$Dn, p=gof_p$ks),
              CvM=list(W=gof_p$W, p=gof_p$cramer),
              nboot =  attr(gof_p, "nboot"),
              boot_success = attr(gof_p, "boot_success")
             )
        )
}


get_edf_cdf <- function(model){
  fun <- function(x, z, lt){
    if(lt){
      length(z[z<x])
    }else{
      length(z[z<=x])
    }
  }

  if("ds" %in% class(model)){
    cdfvalues <- sort(cdf.ds(model)$fitted)
  }else if("io" %in% class(model) |
           "trial" %in% class(model) |
           "rem" %in% class(model) ){
    cdfvalues <- sort(cdf.ds(model$ds)$fitted)
  }else if("io.fi" %in% class(model)){
    data <- model$data
    data <- data[data$object %in% as.numeric(names(model$fitted)), ]
    n <- length(data$distance)/2
    cdfvalues <- rep(0, n)
    for(i in 1:n){
      newdata <- data[data$object %in% as.numeric(names(model$fitted))[i], ]
      cdfvalues[i] <- predict.io.fi(model, newdata=newdata, integrate=TRUE,
                                    int.range=newdata$distance[1])$fitted
      if(model$meta.data$left!=0){
        widthvalue <- predict.io.fi(model, newdata=newdata, integrate=TRUE,
                                    int.range=model$meta.data$width)$fitted
        cdfvalues[i] <- cdfvalues[i]/widthvalue
      }
    }
    if(model$meta.data$left==0){
      cdfvalues <- cdfvalues/predict.io.fi(model, integrate=TRUE)$fitted
    }
    cdfvalues <- sort(cdfvalues)
  }else if("trial.fi" %in% class(model)){
    data <- model$data
    data <- data[data$observer==1&data$object %in%
                   as.numeric(names(model$fitted)), ]
    n <- length(data$distance)
    cdfvalues <- rep(0, n)
    for(i in 1:n){
      newdata <- data[data$object %in% as.numeric(names(model$fitted))[i], ]
      cdfvalues[i] <- predict.trial.fi(model, newdata=newdata, integrate=TRUE,
                                       int.range=newdata$distance[1])$fitted

      if(model$meta.data$left != 0){
        widthvalue <- predict.trial.fi(model, newdata=newdata, integrate=TRUE,
                                       int.range=model$meta.data$width)$fitted
        cdfvalues[i] <- cdfvalues[i]/widthvalue
      }
    }
    if(model$meta.data$left==0){
      cdfvalues <- cdfvalues/predict.trial.fi(model, newdata=data,
                                              integrate=TRUE)$fitted
    }
    cdfvalues <- sort(cdfvalues)
  }else if("rem.fi" %in% class(model)){
    data <- model$data
    data <- data[data$object %in% as.numeric(names(model$fitted)),]
    n <- length(data$distance)/2
    cdfvalues <- rep(0, n)
    for(i in 1:n){
      newdata <- data[data$object %in% as.numeric(names(model$fitted))[i], ]
      cdfvalues[i] <- predict.rem.fi(model, newdata=newdata, integrate=TRUE,
                                     int.range=newdata$distance[1])$fitted

      if(model$meta.data$left!=0){
        widthvalue <- predict.rem.fi(model, newdata=newdata, integrate=TRUE,
                                  int.range=model$meta.data$width)$fitted
        cdfvalues[i] <- cdfvalues[i]/widthvalue
      }
    }
    if(model$meta.data$left==0){
      cdfvalues <- cdfvalues/predict.rem.fi(model, newdata=data,
                                            integrate=TRUE)$fitted
    }
    cdfvalues <- sort(cdfvalues)
  }

  n <- length(cdfvalues)
  lower.edf <- (unlist(sapply(cdfvalues, fun, z=cdfvalues, lt=TRUE)))/n
  upper.edf <- (unlist(sapply(cdfvalues, fun, z=cdfvalues, lt=FALSE)))/n

  return(list(lower.edf = lower.edf,
              upper.edf = upper.edf,
              cdfvalues = cdfvalues,
              n         = n))

}



