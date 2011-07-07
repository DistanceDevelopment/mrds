\name{survey.region.dht}
\alias{survey.region.dht}
\title{Extrapolate Horvitz-Thompson abundance estimates to entire surveyed region...}
\usage{survey.region.dht(Nhat.by.sample, samples, width)}

\description{
  Extrapolate Horvitz-Thompson abundance estimates to
  entire surveyed region :
}
\value{Revised Nhat.by.sample dataframe containing estimates extrapolated
to survey region}
\note{Internal function called by \code{\link{dht}} and related functions.}
\author{Jeff Laake}
\keyword{utility}
\arguments{
  \item{Nhat.by.sample}{dataframe of abundance by sample}
  \item{samples}{samples table}
  \item{width}{transect width}
}
