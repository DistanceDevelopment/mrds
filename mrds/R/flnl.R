#' Log-likelihood computation for distance sampling data
#' 
#' For a specific set of parameter values, it computes and returns the negative
#' log-likelihood for the distance sampling likelihood for distances that are
#' unbinned, binned and a mixture of both.  The function \code{flnl} is the
#' function minimized using \code{\link{optim}} from within
#' \code{\link{ddf.ds}}.  
#' 
#' Most of the computation is in \code{flt.lnl} for line
#' transect data and \code{fpt.lnl} for point count data in which the negative
#' log-likelihood is computed for each observation. \code{flnl} is a wrapper
#' that optionally outputs intermediate results and sums the individual
#' log-likelihood values.
#' 
#' \code{flnl} is the main routine that manipulates the parameters using
#' \code{\link{getpar}} to handle fitting of key, adjustment or all of the
#' parameters.  It then calls \code{flt.lnl} for line transect data or
#' \code{fpt.lnl} for point count data to do the actual computation of the
#' likelihood. Both of those functions can handle a mixture of binned
#' (interval) or unbinned distance measurements.  The probability density
#' function for point counts is \code{fr} and for line transects is \code{fx}.
#' \code{fx}=g(x)/mu (where g(x) is the detection function); whereas,
#' f(r)=r*g(r)/mu where mu in both cases is the normalizing constant.  Both
#' functions are in source code file for \code{link{detfct}} The integral
#' calculations are made with \code{\link{integratedetfct}} and
#' \code{\link{tablecgf}} is used for some detection functions to create a
#' "lookup" table of sorts for standardized integral values that can be scaled
#' much like standard normal distribution.
#' 
#' @aliases flnl flt.lnl fpt.lnl
#' @param fpar parameter values for detection function at which log-likelihood
#'   should be evaluated
#' @param ddfobj distance sampling object
#' @param TCI TRUE if point independence assumed (only relevant for double
#'   observer survey analysis)
#' @param misc.options width-transect width (W); int.range-integration range
#'   for observations; showit-if TRUE shows values of parameters and
#'   log-likelihood; doeachint-if TRUE doesn't use cgftab and does each
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
flnl <- function(fpar, ddfobj, TCI, misc.options, fitting="all")
#
# flt - computes objective function for line transect fitting of grouped/ungrouped distances
#
# Arguments:
#
# fpar            - det fct parameter values
# ddfobj           - distance sampling object
# TCI             - TRUE if point independence
# misc.options
#	width           - transect width
#   int.range       - integration range for observations
#   showit          - if TRUE shows iteration values
#   doeachint       - if TRUE doesn't use cgftab and does each integral instead (for testing only)
#   integral.numeric- if TRUE integral is computed numerically rather than analytically
# fitting         - "key" if only fitting key fct parameters, "adjust" if fitting adjustment function parameters
#                     or "all" to fit both.
#
#  Value:
#
#  lnl - sum of negative log-likelihood values at values of fpar
#
#  Functions Used:  flt.lnl, fpt.lnl, getpar
#
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
      save.shape<-ddfobj$shape$parameters
      ddfobj$shape$parameters=rep(NA,length(ddfobj$shape$parameters))
      pars=getpar(ddfobj)
      fpar[which(is.na(pars))]<-save.shape
    }

    if(!is.null(ddfobj$scale)){
      save.scale<-ddfobj$scale$parameters
      ddfobj$scale$parameters=rep(NA,length(ddfobj$scale$parameters))

      pars=getpar(ddfobj)
      fpar[which(is.na(pars))]<-save.scale
    }
  } 		

  #pars=getpar(ddfobj)
  #pars[which(is.na(pars))]<-fpar[which(is.na(pars))]
  #fpar<-pars

  if(misc.options$point){
    lnl<-sum(fpt.lnl(fpar, ddfobj, TCI, misc.options))
  }else{
    lnl<-sum(flt.lnl(fpar, ddfobj, TCI, misc.options))
  }

  # If iteration results are printed, output log-likelihood
  if(misc.options$showit==3)
    cat("lt lnl = ", lnl,   "\n")    
  
  return(lnl)
}
