\name{gstdint}
\alias{gstdint}
\title{Integral of standardized (scale=1) detection function}
\usage{
  gstdint(x, xmin = NULL, ddfobj, index = NULL,
    select = NULL, width, standardize = TRUE,
    point = FALSE)
}
\arguments{
  \item{x}{vector of lower or upper values for integral
  bound (depending on value of xmin)}

  \item{xmin}{lower bound if specified}

  \item{ddfobj}{distance detection function specification}

  \item{index}{specific data row index}

  \item{select}{logical vector for selection of data
  values}

  \item{width}{truncation width}

  \item{standardize}{logical used to decide whether to
  divide through by the function evaluated at 0}

  \item{point}{logical to determine if point count(TRUE) or
  line transect(FALSE)}
}
\value{
  vector of integral values of detection function with
  scale=1
}
\description{
  Computes the integral of a standardized (scale=1)
  detection function for a specified range
}
\details{
  If \code{xmin} is specified it computes the integral from
  \code{xmin} to each value in \code{x}; otherwise, it
  computes the integral from 0 to \code{x}.
}
\note{
  This is an internal function that is not intended to be
  invoked directly.
}
\author{
  Jeff Laake
}
\seealso{
  \code{\link{tablecgf}}
}
\keyword{utility}

