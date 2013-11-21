\name{pdot.dsr.integrate.logistic}
\alias{pdot.dsr.integrate.logistic}
\title{Compute probability that a object was detected by at least one observer}
\usage{
pdot.dsr.integrate.logistic(right, width, beta, x, integral.numeric, BT, models,
  GAM = FALSE, rem = FALSE, point = FALSE)
}
\arguments{
  \item{right}{???}

  \item{width}{transect width}

  \item{beta}{parameters of logistic detection function}

  \item{x}{data matrix}

  \item{integral.numeric}{???}

  \item{BT}{???}

  \item{models}{list of models including \code{g0model}}

  \item{GAM}{???=FALSE}

  \item{rem}{???=FALSE}

  \item{point}{\code{TRUE} for point transects}
}
\description{
Computes probability that a object was detected by at least
one observer (\code{pdot} or p_.) for a logistic detection
function that contains distance.
}
\author{
Jeff Laake
}

