\name{create.ddfobj}
\alias{create.ddfobj}
\title{Create detection function object}
\usage{
create.ddfobj(model, truncation, xmat, meta.data, initial)
}
\arguments{
  \item{model}{model list with key function and possibly
  adjustment functions, scale formula, and shape formula}

  \item{truncation}{either: a single number giving the
  right truncation for the distances; a 2-vector giving the
  left and right truncations (\code{c(left, right)}); or a
  list with elements \code{left} and \code{right} (at least
  \code{right} must be supplied in list form). Default is
  \code{NULL} which uses the largest observed distance (not
  usually a good idea with unbinned data).}

  \item{xmat}{model data frame}

  \item{meta.data}{list of options}

  \item{initial}{vector of initial values for parameters of
  the detection function}
}
\value{
Distance sampling function object list with elements that
all can be null except type: \item{type}{type of detection
function hn,hr,gamma,unif,logistic} \item{xmat}{model data
frame} \item{intercept.only}{TRUE if scale = ~1 and any
shape formula =~1} \item{cgftab}{table of standardized
integral values} \item{scale}{sublist with elements (can be
NULL i.e., unif key): formula, parameters, design matrix
(dm)} \item{shape}{sublist with elements (power of hazard
rate or gamma) (can be NULL i.e., unif or hn key):formula,
parameters, design matrix (dm)} \item{adjustment}{sublist
with elements (is NULL if no adjustments used):
series,order,scale,parameters} \item{g0}{sublist with
elements (not used at present):formula,parameters, design
matrix(dm), link}
}
\description{
Creates and populates a specific list structure to define a
detection function object and its data. The \code{ddfobj}
is used throughout the package as a calling argument to
various functions.
}
\note{
Internal function not meant to be called by user
}
\author{
Jeff Laake
}
\seealso{
\code{\link{detfct}}, \code{\link{ddf}}
}

