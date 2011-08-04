#' Collection of functions for logistic detection functions
#' 
#' These functions are used to test whether a logistic detection function is a
#' linear function of distance (\code{is.linear.logistic}) or is constant
#' (varies by distance but no other covariates) \code{is.logistic.constant}).
#' Based on these tests, the most appropriate manner for integrating the
#' detection function with respect to distance is chosen.  The integrals are
#' needed to estimate the average detection probability for a given set of
#' covariates.
#' 
#' If the logit is linear in distance then the integral can be computed
#' analytically. If the logit is constant or only varies by distance then only
#' one integral needs to be computed rather than an integral for each
#' observation.
#' 
#' integratelogistic.analytic: as long as distance is a linear function in the
#' logistic the integral can be computed analytically. For reference see
#' integral 526 in CRC Std Math Table 24th ed.
#' 
#' @usage is.linear.logistic(xmat,g0model,zdim,width)
#'        is.logistic.constant(xmat, g0model, width)
#'        integratelogistic.analytic(x, models, beta, width)
#' @aliases logisticdetfct logisticbyx logisticbyz logisticdupbyx
#'   is.linear.logistic is.logistic.constant integratedetfct.logistic
#'   integratelogistic integratelogistic.analytic integratelogisticdup
#'   pdot.dsr.integrate.logistic
#' @param xmat data matrix
#' @param x data matrix (same as xmat)
#' @param g0model logit model
#' @param models list of models including \code{g0model}
#' @param beta parameters of logistic detection function
#' @param zdim number of columns in design matrix
#' @param width transect width
#' @return Logical TRUE if condition holds and FALSE otherwise
#' @author Jeff Laake
#' @keywords utility
is.linear.logistic <-
function(xmat,g0model,zdim,width)
#
#  is.linear.logistic - determines whether the specified logit model is linear for distance
#  If it is linear then the integral can be computed analytically.
#
# Arguments:
#
# xmat     - data
# g0model  - logit model
# zdim     - # of columns in design matrix
# width    - transect width
#
# Value: logical value integral.numeric
#
{
    xmat$distance <- rep(width/2, dim(xmat)[1])
    beta <- rep(1,zdim)
    logit1 <- mean(beta %*% t(setcov(xmat, g0model)$cov))
    xmat$distance <- rep(width, dim(xmat)[1])
    logit2 <- mean(beta %*% t(setcov(xmat, g0model)$cov))
    xmat$distance <- rep(0, dim(xmat)[1])
    logit0 <- mean( beta %*% t(setcov(xmat, g0model)$cov))
        if(logit1-logit0==0)
           integral.numeric=FALSE
        else
       if((logit2-logit0)/(logit1-logit0) <= 2.00001 & (logit2-logit0)/(logit1-logit0) >= 1.99999)
            integral.numeric <- FALSE
       else
            integral.numeric <- TRUE
}
