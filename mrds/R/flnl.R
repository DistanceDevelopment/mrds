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
#' parameters.  It then calls \code{flpt.lnl} to do the actual computation of the
#' likelihood.  The probability density
#' function for point counts is \code{fr} and for line transects is \code{fx}.
#' \code{fx}=g(x)/mu (where g(x) is the detection function); whereas,
#' f(r)=r*g(r)/mu where mu in both cases is the normalizing constant.  Both
#' functions are in source code file for \code{link{detfct}} and are called from
#' distpdf and the integral calculations are made with \code{\link{integratepdf}} and
#' \code{\link{tablecgf}} is used to create a
#' "lookup" table of sorts for standardized integral values that can be scaled
#' much like standard normal distribution.
#' 
#' @aliases flnl flpt.lnl 
#' @param fpar parameter values for detection function at which log-likelihood
#'   should be evaluated
#' @param ddfobj distance sampling object
#' @param misc.options width-transect width (W); int.range-integration range
#'   for observations; showit- 0 to 3 controls level of iteration output; 
#'   doeachint-if TRUE doesn't use cgftab and does each
#'   integral; integral.numeric-if TRUE integral is computed numerically rather
#'   than analytically
#' @param fitting "key" if only fitting key fct parameters, "adjust" if fitting
#'   adjustment function parameters or "all" to fit both
#' @return negative log-likelihood value at the parameter values specified in
#'   \code{fpar}
#' @note These are internal functions used by \code{\link{ddf.ds}} to fit
#'   distance sampling detection functions.  It is not intended for the user to
#'   invoke these functions but they are documented here for completeness.
#' @author Jeff Laake
#' @seealso \code{\link{flt.var}}, \code{\link{detfct}}
#' @keywords utility
flnl <- function(fpar, ddfobj, misc.options, fitting="all")
{
#
#   If iteration results are printed, output parameter values
#   17-Aug-05 jll changed value of showit here and below
#   dlm 05-June-2006 Added an extra level of showit.
  if(misc.options$showit==3)
    cat("\npar = ", fpar,"\n")

  # dlm 27-May-2006 (or some time around then...)
  # During the optimisation we want to make sure that we are keeping the right things
  # constant, so lets do that...

  # Note we do this backwards so that optim doesn't get confused.
  if(fitting=="key"){
    if(!is.null(ddfobj$adjustment)){
      save.adj<-ddfobj$adjustment$parameters
      ddfobj$adjustment$parameters=rep(NA,length(ddfobj$adjustment$parameters))
      pars=getpar(ddfobj)
      fpar[which(is.na(pars))]<-save.adj
    }
  }else if(fitting=="adjust"){
    if(!is.null(ddfobj$shape)){
      save.shape <- ddfobj$shape$parameters
      ddfobj$shape$parameters <- rep(NA,length(ddfobj$shape$parameters))
      pars <- getpar(ddfobj)
      fpar[which(is.na(pars))] <- save.shape
    }

    if(!is.null(ddfobj$scale)){
      save.scale<-ddfobj$scale$parameters
      ddfobj$scale$parameters=rep(NA,length(ddfobj$scale$parameters))

      pars=getpar(ddfobj)
      fpar[which(is.na(pars))]<-save.scale
    }
  } 		
#  compute total negative log-likelihood
   lnl <- sum(flpt.lnl(fpar, ddfobj, misc.options))

  # If iteration results are printed, output log-likelihood
  if(misc.options$showit==3){
    cat("\npar = ", fpar,"\n")
    cat("lt lnl = ", lnl,   "\n")    
  }
  return(lnl)
}
