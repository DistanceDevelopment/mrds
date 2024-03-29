#' Log-likelihood computation for distance sampling data
#'
#' For a specific set of parameter values, it computes and returns the negative
#' log-likelihood for the distance sampling likelihood for distances that are
#' unbinned, binned and a mixture of both.  The function \code{flnl} is the
#' function minimized using \code{\link{optim}} from within
#' \code{\link{ddf.ds}}.
#'
#' Most of the computation is in \code{flpt.lnl} in which the negative
#' log-likelihood is computed for each observation. \code{flnl} is a wrapper
#' that optionally outputs intermediate results and sums the individual
#' log-likelihood values.
#'
#' \code{flnl} is the main routine that manipulates the parameters using
#' \code{\link{getpar}} to handle fitting of key, adjustment or all of the
#' parameters.  It then calls \code{flpt.lnl} to do the actual computation of
#' the likelihood.  The probability density function for point counts is
#' \code{fr} and for line transects is \code{fx}.
#' \code{fx}=g(x)/mu (where g(x) is the detection function); whereas,
#' f(r)=r*g(r)/mu where mu in both cases is the normalizing constant.  Both
#' functions are in source code file for \code{link{detfct}} and are called from
#' \code{distpdf} and the integral calculations are made with
#' \code{\link{integratepdf}}.
#'
#' @aliases flnl flpt.lnl
#' @param fpar parameter values for detection function at which negative
#' log-likelihood should be evaluated
#' @param ddfobj distance sampling object
#' @param misc.options a \code{list} with the following elements: \code{width}
#' transect width; \code{int.range} the integration range for observations;
#' \code{showit} 0 to 3 controls level debug output; \code{integral.numeric} if
#' \code{TRUE} integral is computed numerically rather than analytically;
#' \code{point} is this a point transect?
#' @param fitting character \code{"key"} if only fitting key function
#' parameters, \code{"adjust"} if fitting adjustment parameters or \code{"all"}
#' to fit both
#' @return negative log-likelihood value at the parameter values specified in
#' \code{fpar}
#' @note These are internal functions used by \code{\link{ddf.ds}} to fit
#' distance sampling detection functions.  It is not intended for the user to
#' invoke these functions but they are documented here for completeness.
#' @author Jeff Laake, David L Miller
#' @seealso \code{\link{flt.var}}, \code{\link{detfct}}
#' @keywords utility
flnl <- function(fpar, ddfobj, misc.options, fitting="all"){

  # if we're fitting by key or adjustments only, add in the other parameters
  # to ensure that we can evaluate the likelihood. The other values
  # are in ddfobj (see corresponding code in detfct.fit.opt)
  if(fitting=="key"){
    allvals <- getpar(ddfobj)
    parind <- getpar(ddfobj, index=TRUE)
    allvals[1:parind[2]] <- fpar
    fpar <- allvals
  }else if(fitting=="adjust"){
    allvals <- getpar(ddfobj)
    parind <- getpar(ddfobj, index=TRUE)
    allvals[(parind[2]+1):parind[3]] <- fpar
    fpar <- allvals
  }

  #  compute total negative log-likelihood
  nll <- sum(flpt.lnl(fpar, ddfobj, misc.options))

  # If iteration results are printed, output
  # log-likelihood and parameter values
  if(misc.options$showit >= 3){
    cat("par = ", fpar,"\n")
    cat("nll = ", nll,   "\n")
  }
  return(nll)
}
