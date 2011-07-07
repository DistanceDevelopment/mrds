\name{test.breaks}
\alias{test.breaks}
\title{Test validity for histogram breaks(cutpoints)...}
\usage{test.breaks(breaks, left, width)}

\description{
  Test validity for histogram breaks(cutpoints) Determines
  whether user specified breaks for histograms are properly
  ordered and match the left and right truncation.
}
\value{vector of breaks modified to be valid if necessary}
\author{Jeff Laake}
\arguments{
  \item{breaks}{vector of cutpoints (breaks) for distance histogram}
  \item{left}{left truncation value}
  \item{width}{right truncation value; either radius of point count or
half-width of transect}
}
