\name{ddf.gof}
\alias{ddf.gof}
\alias{gof.ds}
\alias{gof.io}
\alias{gof.io.fi}
\alias{gof.rem}
\alias{gof.rem.fi}
\alias{gof.trial}
\alias{gof.trial.fi}
\title{Goodness of fit tests for distance sampling models}
\usage{
  ddf.gof(model, breaks = NULL, nc = NULL, qq = TRUE, ...)
}
\arguments{
  \item{model}{ddf model object}

  \item{breaks}{Cutpoints to use for binning data}

  \item{nc}{Number of distance classes}

  \item{qq}{Flag to indicate whether quantile-quantile plot
  is desired}

  \item{\dots}{Graphics parameters to pass into qqplot
  function}
}
\value{
  List of class 'ddf.gof' containing \item{chi-square
  }{Goodness of fit test statistic} \item{df }{Degrees of
  freedom associated with test statistic} \item{p-value
  }{Significance level of test statistic}
}
\description{
  Generic function that computes chi-square goodness of fit
  test for ddf models
}
\author{
  Jeff Laake
}
\seealso{
  \code{\link{qqplot.ddf}}
}
\keyword{utility}

