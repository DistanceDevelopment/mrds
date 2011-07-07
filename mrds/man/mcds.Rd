\name{mcds}
\alias{mcds}
\title{MCDS function definition...}
\usage{mcds(formula, key, adj.series, adj.order=c(NULL),
    adj.scale="width", shape.formula=~1)}

\description{
  MCDS function definition Creates model formula list for
  multiple covariate distance sampling using values
  supplied in call to \code{\link{ddf}}
}
\value{A formula list used to define the detection function model
\item{fct}{string "mcds"} \item{key}{key function string}
\item{adj.series}{adjustment function string} \item{adj.order}{adjustment
function orders} \item{adj.scale}{adjustment function scale type}
\item{formula}{formula for scale function} \item{shape.formula}{formula
for shape function}}
\author{Jeff Laake; Dave Miller}
\keyword{utility}
\arguments{
  \item{formula}{formula for scale function}
  \item{key}{string identifying key function (currently either "hn"
(half-normal),"hr" (hazard-rate), "unif" (uniform) or "gamma" (gamma
distribution)}
  \item{adj.series}{string identifying adjustment functions cos (Cosine), herm
(Hermite polynomials), poly (simple polynomials) or NULL}
  \item{adj.order}{vector of order of adjustment terms to include}
  \item{adj.scale}{whether to scale the adjustment terms by "width" or "scale"}
  \item{shape.formula}{formula for shape function}
}
