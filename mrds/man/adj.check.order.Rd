\name{adj.check.order}
\alias{adj.check.order}
\title{Check order of adjustment terms}
\usage{
  adj.check.order(adj.series, adj.order)
}
\arguments{
  \item{adj.series}{Adjustment series used
  ('\code{cos}','\code{herm}','\code{poly}')}

  \item{adj.order}{Integer to check}
}
\value{
  Logical value, true if there have been errors, false
  otherwise.
}
\description{
  'adj.check.order' checks that the Cosine, Hermite or
  simple polynomials are of the correct order.
}
\details{
  Only even functions are allowed as adjustment terms. Also
  Hermite polynomials must be of degree at least 4 and
  Cosine of order at least 3. If incorrect terms are
  supplied then \code{errors} is called.
}
\author{
  David Miller
}
\references{
  S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake.
  1993.  Robust Models. In: Distance Sampling, eds.
  S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake.
  Chapman & Hall.
}
\seealso{
  \code{\link{adjfct.cos}}, \code{\link{adjfct.poly}},
  \code{\link{adjfct.herm}}, \code{\link{detfct}},
  \code{\link{mcds}}, \code{\link{cds}}
}
\keyword{methods}

