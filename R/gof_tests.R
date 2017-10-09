#' Goodness of fit testing for detection functions
#'
#' Formal goodness of fit testing for detection function models using Kolmogorov-Smirnov and Cramer-von Mises tests. Both tests are based on looking at the quantile-quantile plot produced by \code{\link{qqplot.ddf}} and deviations from the line x=y.
#'
#' The Kolmogorov-Smirnov test asks the question “what’s the largest vertical distance between a point and the y = x line?” It uses this distance as a statistic to test the null hypothesis that the samples (EDF and CDF in our case) are from the same distribution (and hence our model fits well). If the deviation between the y = x line and the points is too large we reject the null hypothesis and say the model doesn’t have a good fit.
#'
#' Rather than looking at the single biggest difference between the y = x line and the points in the Q-Q plot, we might prefer to think about all the differences between line and points, since there may be many smaller differences that we want to take into account rather than looking for one large deviation. Its null hypothesis is the same, but the statistic it uses is the sum of the deviations from each of the point to the line.
#
#' @section Details:

#' Note that a bootstrap procedure is required to ensure that the p-values from the procedure are correct as the we are comparing the cumulative distribution function (CDF) and empirical distribution function (EDF) and we have estimated the parameters of the detection function.
#'
#' @param model a fitted model object
#' @param nboot number of bootstraps to do
#'
#' @return list of p-values for the two tests (\code{ks.p}, \code{cramer.p}) an dtest statistics (\code{Dn} for K-S and \code{W} for C-vM).
#' @author David L Miller
#'
gof_tests <- function(model, nboot=100){

  # storage
  ks.boot <- rep(0, nboot)
  cramer.boot <- rep(0, nboot)

  for (i in 1:nboot) {
    # simulate data from the model
    # fit a model
    refit <- sample_ddf(model)

    if(!inherits(refit, "try-error")){

      ## calculate test statistics and store them
      edf_cdf <- get_edf_cdf(model)
      # Kolmogorov-Smirnov
      ks.boot[i] <- max(c(abs(edf_cdf$lower.edf - edf_cdf$cdfvalues),
                          abs(edf_cdf$upper.edf - edf_cdf$cdfvalues)))
      # Cramer-von Mises
      cramer.boot[i] <- 1/(12*edf_cdf$n) + sum((edf_cdf$cdfvalues -
                                                ((1:edf_cdf$n)-.5)/edf_cdf$n)^2)

    }
  }

  # now calculate for the model
  Dn <- max(c(abs(edf_cdf$lower.edf - edf_cdf$cdfvalues),
              abs(edf_cdf$upper.edf - edf_cdf$cdfvalues)))
  W <- 1/(12*edf_cdf$n) + sum((edf_cdf$cdfvalues - ((1:edf_cdf$n)-.5)/edf_cdf$n)^2)

  # do the boostrap test
  ks.p <- mean(Dn <= ks.boot, na.rm=TRUE)
  cramer.p <- mean(W <= cramer.boot, na.rm=TRUE)

  return(list(ks=ks.p, cramer=cramer.p,
              Dn=Dn, W=W))
}
