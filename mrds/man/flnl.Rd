\name{flnl}
\alias{flnl}
\alias{flpt.lnl}
\title{Log-likelihood computation for distance sampling data}
\usage{
  flnl(fpar, ddfobj, misc.options, fitting = "all")
}
\arguments{
  \item{fpar}{parameter values for detection function at
  which log-likelihood should be evaluated}

  \item{ddfobj}{distance sampling object}

  \item{misc.options}{width-transect width (W);
  int.range-integration range for observations; showit- 0
  to 3 controls level of iteration output; doeachint-if
  TRUE doesn't use cgftab and does each integral;
  integral.numeric-if TRUE integral is computed numerically
  rather than analytically}

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
  Most of the computation is in \code{flpt.lnl} in which
  the negative log-likelihood is computed for each
  observation. \code{flnl} is a wrapper that optionally
  outputs intermediate results and sums the individual
  log-likelihood values.

  \code{flnl} is the main routine that manipulates the
  parameters using \code{\link{getpar}} to handle fitting
  of key, adjustment or all of the parameters.  It then
  calls \code{flpt.lnl} to do the actual computation of the
  likelihood.  The probability density function for point
  counts is \code{fr} and for line transects is \code{fx}.
  \code{fx}=g(x)/mu (where g(x) is the detection function);
  whereas, f(r)=r*g(r)/mu where mu in both cases is the
  normalizing constant.  Both functions are in source code
  file for \code{link{detfct}} and are called from distpdf
  and the integral calculations are made with
  \code{\link{integratepdf}} and \code{\link{tablecgf}} is
  used to create a "lookup" table of sorts for standardized
  integral values that can be scaled much like standard
  normal distribution.
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

