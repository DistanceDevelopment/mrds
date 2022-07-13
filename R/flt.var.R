#' Hessian computation for fitted distance detection function model parameters
#'
#' Computes hessian to be used for variance-covariance matrix.  The hessian is
#' the outer product of the vector of first partials (see pg 62 of Buckland et
#' al 2002).
#'
#' @param ddfobj distance sampling object
#' @param misc.options width-transect width (W); int.range-integration range
#'   for observations; showit-0 to 3 controls level of iteration printing;
#'   integral.numeric-if TRUE integral is computed numerically rather
#'   than analytically
#' @return variance-covariance matrix of parameters in the detection function
#' @note This is an internal function used by \code{\link{ddf.ds}} to fit
#'   distance sampling detection functions.  It is not intended for the user to
#'   invoke this function but it is documented here for completeness.
#' @author Jeff Laake and David L Miller
#' @seealso \code{\link{flnl}},\code{\link{flpt.lnl}},\code{\link{ddf.ds}}
#' @references Buckland et al. 2002
#' @keywords utility
flt.var <- function(ddfobj, misc.options){

  fpar1 <- getpar(ddfobj)
  #   Compute first partial (numerically) of log(f(y)) for each observation
  #   for each parameter and store in parmat (n by length(fpar))
  nflpt.lnl <- function(...) -flpt.lnl(...)
  parmat <- numDeriv::jacobian(nflpt.lnl, x=fpar1, ddfobj=ddfobj,
                               misc.options=misc.options)

  # Compute varmat using first partial approach (pg 62 of Buckland et al 2002)
  varmat <- t(parmat) %*% parmat

  # could get the "true" (2nd deriv) hessian using the below
  #nflpt.lnl <- function(...) sum(flpt.lnl(...))
  #varmat <- numDeriv::hessian(nflpt.lnl, x=fpar1, ddfobj=ddfobj,
  #                             misc.options=misc.options)

  return(varmat)
}
