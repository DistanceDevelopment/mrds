\name{flt.var}
\alias{flt.var}
\title{Hessian computation for fitted distance detection function model parameters}
\usage{
  flt.var(ddfobj, TCI, misc.options)
}
\arguments{
  \item{ddfobj}{distance sampling object}

  \item{TCI}{TRUE if point independence assumed (only
  relevant for double observer survey analysis)}

  \item{misc.options}{width-transect width (W);
  int.range-integration range for observations; showit-if
  TRUE shows values of parameters and log-likelihood;
  doeachint-if TRUE doesn't use cgftab and does each
  integral; integral.numeric-if TRUE integral is computed
  numerically rather than analytically}
}
\value{
  variance-covariance matrix of parameters in the detection
  function
}
\description{
  Computes hessian to be used for variance-covariance
  matrix.  The hessian is the outer product of the vector
  of first partials (see pg 62 of Buckland et al 2002).
}
\note{
  This is an internal function used by \code{\link{ddf.ds}}
  to fit distance sampling detection functions.  It is not
  intended for the user to invoke this function but it is
  documented here for completeness.
}
\author{
  Jeff Laake
}
\references{
  Buckland et al. 2002
}
\seealso{
  \code{\link{flnl}},\code{\link{flt.lnl}},\code{\link{ddf.ds}}
}
\keyword{utility}

