\name{setinitial.ds}
\title{Set initial values for detection function based on distance sampling...}
\usage{setinitial.ds(ddfobj,width,initial,point)
sethazard(ddfobj,dmat,width)}

\description{
  Set initial values for detection function based on
  distance sampling For a given detection function, it
  computes the initial values for the parameters including
  scale and shape parameters and adjustment function
  parameters if any.  If there are user-defined initial
  values only the parameters not specified by the user are
  computed.
}
\alias{setinitial.ds}
\alias{sethazard}
\value{\item{scale}{vector of initial scale parameter values}
\item{shape}{vector of initial shape parameter values}
\item{adjustment}{vector of initial adjustment function parameter values}}
\author{Jeff Laake, Dave Miller}
\arguments{
  \item{ddfobj}{distance detection function object}
  \item{width}{half-width of transect or radius of point count}
  \item{initial}{list of user-defined initial values with possible elements
scale,shape,adjustment}
  \item{point}{if TRUE, point count data; otherwise, line transect data}
  \item{dmat}{xmat from ddfobj}
}
