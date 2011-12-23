\name{qqplot.ddf}
\alias{pcramer}
\alias{pks}
\alias{qqplot.ddf}
\title{Q-Q plot, KS and CVM goodness of fit tests for distance detection functions}
\usage{
  qqplot.ddf(model,...)

  pcramer(q, eps = 1e-05)

  pks(Dn,n)
}
\arguments{
  \item{model}{fitted distance detection function model
  object}

  \item{n}{sample size}

  \item{Dn}{K-S statistic}

  \item{q}{CvM statistic}

  \item{eps}{small value that controls accuracy of p-value
  computation}

  \item{\dots}{unspecified arguments passed to plot}
}
\value{
  A list of goodness of fit related values:
  \item{edf}{matrix of lower and upper empirical
  distribution function values} \item{cdf}{fitted
  cumulative distribution function values} \item{ks}{list
  with K-S statistic (\code{Dn}) and p-value (\code{p})}
  \item{CvM}{list with CvM statistic (\code{W}) and p-value
  (\code{p})}
}
\description{
  Constructs a Q-Q plot for fitted model as a graphical
  picture of goodness of fit and computes K-S and
  Cramer-VonMises goodness of fit tests for distance
  sampling models based on single observer survey and
  double observer survey with independent observer (io) and
  trial configurations.
}
\details{
  \code{pks} computes the p-value for the K-S test.  The
  function \code{pcramer} was taken from the coda package.
  It computes the p-value for the CvM test.  Both
  \code{pks} and \code{pcramer} are used in
  \code{qqplot.ddf} and need not be called by user.
  \code{qqplot.ddf} is called from \code{ddf.gof} to
  evaluate model goodness of fit.
}
\author{
  Jeff Laake
}
\references{
  Burnham, K.P., S.T. Buckland, J.L. Laake, D.L. Borchers,
  T.A.  Marques, J.R.B. Bishop, and L. Thomas. 2004.
  Further topics in distance sampling. pp: 385-389. In:
  Advanced Distance Sampling, eds. S.T. Buckland,
  D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers,
  and L. Thomas.  Oxford University Press.
}
\seealso{
  \code{\link{ddf.gof}}, \code{\link{cdf.ds}}
}
\keyword{utility}

