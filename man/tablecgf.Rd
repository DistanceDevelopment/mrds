\name{tablecgf}
\alias{tablecgf}
\title{Spline approximation to scale dependent integration of detection function}
\usage{
tablecgf(ddfobj, width, standardize = TRUE, point = FALSE)
}
\arguments{
  \item{ddfobj}{distance detection function specification}

  \item{width}{truncation width}

  \item{standardize}{logical used to decide whether to
  divide through by the function evaluated at 0}

  \item{point}{logical to determine if point count(TRUE) or
  line transect(FALSE)}
}
\value{
a smooth.spline result object
}
\description{
Computes spline approximation to cumulative detection
function for interpolation of scale dependent integrations
of the detection function.
}
\details{
This is an internal function used to speed up integration
of the detection function with scale covariates. The
detection function is integrated at a series of points from
0 to W and then a spline is fitted to the computed values
which are cumulative (integral from 0 to x < integral from
0 to x+dx fr dx>0).  The spline is then used to predict
values of the integral which depend on the scale which can
depend on observation specific covariates.
}
\note{
This is an internal function that is not intended to be
invoked directly.
}
\author{
Jeff Laake
}
\seealso{
\code{\link{gstdint}}
}
\keyword{utility}

