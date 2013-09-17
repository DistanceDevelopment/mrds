\name{calcp.mrds}
\alias{calcp.mrds}
\title{Average detection probability calculation- DEPRECATED FUNCTION NO LONGER USED}
\usage{
  calcp.mrds(dpformula, dplink, dppars, dpdata, vname,
    lower = 0, upper, divisions = 30, type = TRUE,
    objname = "sighting", obsname = "plat")
}
\arguments{
  \item{dpformula}{valid formula for binary regression
  model detection function}

  \item{dplink}{link function}

  \item{dppars}{vector of detection function parameters,
  named in agreement with dpformula, and in same order (not
  sure if order is crucial)}

  \item{dpdata}{data frame in format for ddf() mrds
  analysis}

  \item{vname}{name of distance variable for x-axis (would
  be "distance" in usual ddf() data frame)}

  \item{lower}{lower bound for plotting on x-axis}

  \item{upper}{upper bound for plotting on x-axis}

  \item{divisions}{number of divisions to use for numerical
  integration}

  \item{type}{"line" for line transects, "point" for point
  transects (affects how integration is done)}

  \item{objname}{name of variable denoting
  detection/sighting/object number (would be "object" in
  usual ddf() data frame)}

  \item{obsname}{name of variable denoting
  observer/platform number (1 or 2) (would be "observer" in
  usual ddf() data frame)}
}
\value{
  \item{x}{distances at which probabilities are calculated}
  \item{p1}{probability for observer 1}
  \item{p2}{probability for observer 2}
  \item{p3}{probability seen by both observers}
  \item{p.}{probability seen by at least one observer}
}
\description{
  Detection probability at a distance can vary depending on
  the other ancillary covariates collected for an
  observation.  To plot a single detection line, it is
  necessary to compute the expected probability value at
  distance across the estimated proportions of the
  covariates in the population.  This function does that
  calculation by observer, for duplicates and cummulative
  for both observers.
}
\author{
  Who wrote this?
}
\references{
  See Laake and Borchers (2004) chapter 6 in Advanced
  Distance Sampling by Buckland et al.
}

