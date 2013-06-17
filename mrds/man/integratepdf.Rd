\name{integratepdf}
\alias{integratepdf}
\title{Numerically integrate pdf of observed distances over specified ranges}
\usage{
  integratepdf(ddfobj, select, width, int.range,
    doeachint = FALSE, standardize = TRUE, point = FALSE)
}
\arguments{
  \item{ddfobj}{distance detection function specification}

  \item{select}{logical vector for selection of data
  values}

  \item{width}{truncation width}

  \item{int.range}{integration range matrix; vector is
  converted to matrix}

  \item{doeachint}{logical that specifies whether each
  observation integral should be computed numerically}

  \item{standardize}{logical used to decide whether to
  divide through by the function evaluated at 0}

  \item{point}{logical to determine if point count(TRUE) or
  line transect(FALSE)}
}
\value{
  vector of integral values - one for each observation
}
\description{
  Computes integral of pdf of observed distances over x for
  each observation.  The method of computation depends on
  argument switches set and the type of detection function.
}
\details{
  If either doeachint is set or there is only one integral
  then they are computed using integrate; otherwise, it
  uses the cgftab which is a spline fitted to a table of
  standardized integrals and the value is interpolated from
  the spline for each observation.
}
\author{
  Jeff Laake
}
\keyword{utility}

