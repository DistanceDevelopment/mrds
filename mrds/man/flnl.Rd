\name{flnl}
\alias{flnl}
\alias{flt.lnl}
\alias{fpt.lnl}
\title{Log-likelihood computation for distance sampling data}
\usage{
  flnl(fpar, ddfobj, TCI, misc.options, fitting = "all")
}
\arguments{
  \item{fpar}{parameter values for detection function at
  which log-likelihood should be evaluated}

  \item{ddfobj}{distance sampling object}

  \item{TCI}{TRUE if point independence assumed (only
  relevant for double observer survey analysis)}

  \item{misc.options}{width-transect width (W);
  int.range-integration range for observations; showit-if
  TRUE shows values of parameters and log-likelihood;
  doeachint-if TRUE doesn't use cgftab and does each
  integral; integral.numeric-if TRUE integral is computed
  numerically rather than analytically}

  \item{fitting}{"key" if only fitting key fct parameters,
  "adjust" if fitting adjustment function parameters or
  "all" to fit both}
}
\value{
  negative log-likelihood value at the parameter values
  specified in \code{fpar}
}
\description{
  For a specific set of parameter values, it computes and
  returns the negative log-likelihood for the distance
  sampling likelihood for distances that are unbinned,
  binned and a mixture of both.  The function \code{flnl}
  is the function minimized using \code{\link{optim}} from
  within \code{\link{ddf.ds}}.
}
\details{
  Most of the computation is in \code{flt.lnl} for line
  transect data and \code{fpt.lnl} for point count data in
  which the negative log-likelihood is computed for each
  observation. \code{flnl} is a wrapper that optionally
  outputs intermediate results and sums the individual
  log-likelihood values.

  \code{flnl} is the main routine that manipulates the
  parameters using \code{\link{getpar}} to handle fitting
  of key, adjustment or all of the parameters.  It then
  calls \code{flt.lnl} for line transect data or
  \code{fpt.lnl} for point count data to do the actual
  computation of the likelihood. Both of those functions
  can handle a mixture of binned (interval) or unbinned
  distance measurements.  The probability density function
  for point counts is \code{fr} and for line transects is
  \code{fx}. \code{fx}=g(x)/mu (where g(x) is the detection
  function); whereas, f(r)=r*g(r)/mu where mu in both cases
  is the normalizing constant.  Both functions are in
  source code file for \code{link{detfct}} The integral
  calculations are made with \code{\link{integratedetfct}}
  and \code{\link{tablecgf}} is used for some detection
  functions to create a "lookup" table of sorts for
  standardized integral values that can be scaled much like
  standard normal distribution.
}
\note{
  These are internal functions used by \code{\link{ddf.ds}}
  to fit distance sampling detection functions.  It is not
  intended for the user to invoke these functions but they
  are documented here for completeness.
}
\author{
  Jeff Laake
}
\seealso{
  \code{\link{flt.var}}, \code{\link{detfct}}
}
\keyword{utility}

