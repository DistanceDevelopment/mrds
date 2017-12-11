#' Goodness of fit tests for distance sampling models
#'
#' Generic function that computes chi-square goodness of fit test for detection function models with binned data and Cramer-von Mises and Kolmogorov-Smirnov (if \code{ks=TRUE})tests for exact distance data. By default a Q-Q plot is generated for exact data (and can be suppressed using the \code{qq=FALSE} argument).
#'
#' Formal goodness of fit testing for detection function models using Kolmogorov-Smirnov and Cramer-von Mises tests. Both tests are based on looking at the quantile-quantile plot produced by \code{\link{qqplot.ddf}} and deviations from the line x=y.
#'
#' The Kolmogorov-Smirnov test asks the question "what's the largest vertical distance between a point and the y=x line?" It uses this distance as a statistic to test the null hypothesis that the samples (EDF and CDF in our case) are from the same distribution (and hence our model fits well). If the deviation between the y=x line and the points is too large we reject the null hypothesis and say the model doesn't have a good fit.
#'
#' Rather than looking at the single biggest difference between the y=x line and the points in the Q-Q plot, we might prefer to think about all the differences between line and points, since there may be many smaller differences that we want to take into account rather than looking for one large deviation. Its null hypothesis is the same, but the statistic it uses is the sum of the deviations from each of the point to the line.
#
#' @section Details:
#'
#' Note that a bootstrap procedure is required for the Kolmogorov-Smirnov test to ensure that the p-values from the procedure are correct as the we are comparing the cumulative distribution function (CDF) and empirical distribution function (EDF) and we have estimated the parameters of the detection function. The \code{nboot} parameter controls the number of bootstraps to use. Set to \code{0} to avoid computing bootstraps (much faster but with no Kolmogorov-Smirnov results, of course).
#'
#' @aliases ddf.gof gof.io gof.io.fi gof.trial gof.trial.fi gof.rem gof.rem.fi
#' @export
#' @param model model object
#' @param breaks Cutpoints to use for binning data
#' @param nc Number of distance classes
#' @param qq Flag to indicate whether quantile-quantile plot is desired
#' @param nboot number of replicates to use to calculate p-values for the Kolmogorov-Smirnov goodness of fit test statistics
#' @param ks perform the Kolmogorov-Smirnov test (this involves many bootstraps so can take a while)
#' @param \dots Graphics parameters to pass into qqplot function
#' @return List of class \code{ddf.gof} containing \item{chi-square }{Goodness of fit test statistic} \item{df}{Degrees of freedom associated with test statistic} \item{p-value }{Significance level of test statistic}
#' @author Jeff Laake
#' @seealso \code{\link{qqplot.ddf}}
#' @keywords utility
ddf.gof <- function(model, breaks=NULL, nc=NULL, qq=TRUE, nboot=100, ks=FALSE,
                    ...){
  # Functions Used: gof.ds, gof.io, gof.io.fi, gof.trial,
  #                 gof.trial.fi, qqplot.df

  if(!is.null(breaks)){
    breaks <- test.breaks(breaks,model$meta.data$left,model$meta.data$width)
    nc <- length(breaks)-1
  }else if(model$meta.data$binned){
    breaks <- model$meta.data$breaks
  }

  # call method specific function
  result <- switch(model$method,
                   ds       = gof.ds(model, breaks, nc),
                   io       = gof.io(model, breaks, nc),
                   io.fi    = gof.io.fi(model, breaks, nc),
                   trial    = gof.trial(model, breaks, nc),
                   trial.fi = gof.trial.fi(model, breaks, nc),
                   rem      = gof.rem(model, breaks, nc),
                   rem.fi   = gof.rem.fi(model, breaks, nc))

  if(!model$meta.data$binned){
    result <- list(chisquare=result,
                   dsgof=qqplot.ddf(model, plot=qq, ks=ks, nboot=nboot, ...))
  }else{
    result <- list(chisquare=result)
  }

  class(result) <- c("ddf.gof")
  return(result)
}
