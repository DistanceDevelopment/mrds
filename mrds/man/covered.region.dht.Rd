\name{covered.region.dht}
\alias{covered.region.dht}
\title{Covered region estimate of abundance from Horvitz-Thompson-like estimator...}
\usage{covered.region.dht(obs, samples, group)}

\description{
  Covered region estimate of abundance from
  Horvitz-Thompson-like estimator Computes H-T abundance
  within covered region by sample.
}
\value{Nhat.by.sample - dataframe of abundance by sample}
\note{Internal function called by \code{\link{dht}} and related functions}
\author{Jeff Laake}
\keyword{utility}
\arguments{
  \item{obs}{observations table}
  \item{samples}{samples table}
  \item{group}{if TRUE compute abundance of group otherwise abundance of
individuals}
}
