\name{adj.check.order}
\alias{adj.check.order}
\title{Check order of adjustment terms...}
\usage{adj.check.order(adj.series, adj.order)}

\description{
  Check order of adjustment terms
}

\details{
  'adj.check.order' checks that the Cosine, Hermite or
  simple polynomials are of the correct order.

  Only even functions are allowed as adjustment terms. Also
  Hermite polynomials must be of degree at least 4 and
  Cosine of order at least 3. If incorrect terms are
  supplied then \code{errors} is called.
}
\value{Logical value, true if there have been errors, false otherwise.}
\author{David Miller}
\seealso{\code{\link{adjfct.cos}}, \code{\link{adjfct.poly}},
\code{\link{adjfct.herm}}, \code{\link{detfct}}, \code{\link{mcds}},
\code{\link{cds}}}
\references{S.T.Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake. 1993.
Robust Models. In: Distance Sampling, eds. S.T.Buckland, D.R.Anderson,
K.P. Burnham, J.L. Laake. Chapman & Hall.}
\keyword{methods}
\arguments{
  \item{adj.series}{Adjustment series used
('\code{cos}','\code{herm}','\code{poly}')}
  \item{adj.order}{Integer to check}
}
