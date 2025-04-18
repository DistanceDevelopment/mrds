% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integratepdf.grad.R
\name{integratepdf.grad}
\alias{integratepdf.grad}
\title{Numerically integrates the non-normalised pdf or the detection function of 
observed distances over specified ranges.}
\usage{
integratepdf.grad(
  par.index,
  ddfobj,
  int.range,
  width,
  standardize = FALSE,
  point = FALSE,
  left = 0,
  pdf.based = TRUE
)
}
\arguments{
\item{par.index}{the index of the parameter of interest}

\item{ddfobj}{the ddf object}

\item{int.range}{vector with the lower and upper bound of the integration}

\item{width}{the truncation width}

\item{standardize}{TRUE if the non-standardised detection function should
be integrated. Only implemented for standardize = FALSE, so users should not
touch this argument and it can probably be removed.}

\item{point}{are the data from point transects (TRUE) or line transects 
(FALSE).}

\item{left}{the left truncation. Defaults to zero.}

\item{pdf.based}{evaluate the non-normalised pdf or the detection function? 
Default is TRUE.}
}
\description{
Gradient of the integral of the detection function, i.e., d beta/d theta in 
the documentation. This gradient of the integral is the same as the integral
of the gradient, thanks to Leibniz integral rule.
}
\details{
For internal use only -- not to be called by \code{mrds} or \code{Distance}
users directly.
}
\author{
Felix Petersma
}
