\name{average.line}
\alias{average.line}
\title{Average detection function line for plotting}
\usage{
  average.line(finebr, obs, model)
}
\arguments{
  \item{finebr}{set of fine breaks in distance over which
  detection function values are averaged and plotted}

  \item{obs}{value of observer for averaging (1-2
  individual observers; 3 duplicates; 4 pooled observation
  team)}

  \item{model}{ddf model object}
}
\value{
  list with 2 elements \item{xgrid}{vector of gridded
  distance values} \item{values }{vector of average
  detection fct values at the xgrid values }
}
\description{
  For models with covariates the detection probability for
  each observation can vary.  This function computes an
  average value for a set of distances to plot an average
  line to graphically represent the fitted model in plots
  that compare histograms and the scatter of individual
  estimated detection probabilities.
}
\note{
  Internal function called from plot functions for ddf
  objects
}
\author{
  Jeff Laake
}
\keyword{utility}

