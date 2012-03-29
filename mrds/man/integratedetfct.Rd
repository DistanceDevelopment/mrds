\name{integratedetfct}
\alias{integratedetfct}
\title{Numerically integrate distance detection function over specified range}
\usage{
  integratedetfct(ddfobj, select, width, int.range,
    doeachint = FALSE, standardize = TRUE, point = FALSE)
}
\arguments{
  \item{ddfobj}{distance detection function specification}

  \item{select}{logical vector for selection of data
  values}

  \item{width}{truncation width}

  \item{int.range}{integration range}

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
  Computes integral of detection function over x for each
  observation.  The method of computation depends on
  argument switches set and the type of detection function.
}
\details{
  If there are only intercepts in the scale model then the
  integrals can be computed once and they apply to all
  observations.  If the scale model includes variables
  other than intercepts then an integral must be computed
  for each observation. Unless \code{doeachint} is set to
  TRUE, this is done by creating a table of standardized
  integrals and interpolating the value for each
  observation.  For the half-normal detection function the
  table is computed at the very beginning and does not need
  to be recomputed.  For the hazard rate the table needs to
  be recomputed for each iteration because the shape
  function changes.  For small sample sizes it may be more
  efficient to set doeachint=TRUE. For full independence
  models, integration of p_1(x)*p_2(x) cannot use the table
  approach (except for the half-normal but it has not been
  implemented) so each integral is computed. Except for
  this latter integral, the integral can be computed by
  integrating a standardized detection function and then
  multiplying by the scale for that observation. This is a
  simple substitution approach for integration. y=x/scale.
  However, it does not work for the product because each
  has a separate scale; although a common scale could be
  found algebraically at least for the half-normal (see
  Laake (1999)).
}
\author{
  Jeff Laake
}
\keyword{utility}

