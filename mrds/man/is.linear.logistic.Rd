\name{is.linear.logistic}
\alias{integratedetfct.logistic}
\alias{integratelogistic}
\alias{integratelogistic.analytic}
\alias{integratelogisticdup}
\alias{is.linear.logistic}
\alias{is.logistic.constant}
\alias{logisticbyx}
\alias{logisticbyz}
\alias{logisticdetfct}
\alias{logisticdupbyx}
\alias{pdot.dsr.integrate.logistic}
\title{Collection of functions for logistic detection functions}
\usage{
  is.linear.logistic(xmat,g0model,zdim,width)
    is.logistic.constant(xmat, g0model, width)
    integratelogistic.analytic(x, models, beta, width)
}
\arguments{
  \item{xmat}{data matrix}

  \item{x}{data matrix (same as xmat)}

  \item{g0model}{logit model}

  \item{models}{list of models including \code{g0model}}

  \item{beta}{parameters of logistic detection function}

  \item{zdim}{number of columns in design matrix}

  \item{width}{transect width}
}
\value{
  Logical TRUE if condition holds and FALSE otherwise
}
\description{
  These functions are used to test whether a logistic
  detection function is a linear function of distance
  (\code{is.linear.logistic}) or is constant (varies by
  distance but no other covariates)
  \code{is.logistic.constant}). Based on these tests, the
  most appropriate manner for integrating the detection
  function with respect to distance is chosen.  The
  integrals are needed to estimate the average detection
  probability for a given set of covariates.
}
\details{
  If the logit is linear in distance then the integral can
  be computed analytically. If the logit is constant or
  only varies by distance then only one integral needs to
  be computed rather than an integral for each observation.

  integratelogistic.analytic: as long as distance is a
  linear function in the logistic the integral can be
  computed analytically. For reference see integral 526 in
  CRC Std Math Table 24th ed.
}
\author{
  Jeff Laake
}
\keyword{utility}

