% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rescale_pars.R
\name{rescale_pars}
\alias{rescale_pars}
\title{Calculate the parameter rescaling for parameters associated with covariates}
\usage{
rescale_pars(initialvalues, ddfobj)
}
\arguments{
\item{initialvalues}{starting values for the optimisation}

\item{ddfobj}{detection function object}
}
\description{
This will calculate the rescaling needed when covariates to be included in
the scale of the detection function are "too big". Based on code from
\code{\link[optimx]{optimx}}.
}
\details{
Derivative-free methods like nlminb are sensitive to the parameters being
poorly scaled. This can also cause problems for quasi-Newton methods too (at
least, bad scaling won't _help_ the optimisation). So here we rescale the
parameters if necessary (unless we already got scaling from control)
}
\author{
David L Miller
}
