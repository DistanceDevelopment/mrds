\name{foptim}
\title{Alternate detection function fitting code...}
\usage{foptim(initialvalues, bounds, detfct.options, misc.options)}

\description{
  Alternate detection function fitting code A wrapper
  function for foptim.dll -- the Fortran optimisation
  routine. Most of the code juggles the data before making
  the call to ensure Fortran doesn't fail.
}

\details{
  foptim.dll is a re-implementation of the CDS engine found
  in Distance. It uses the DNCONG/DNCONF routines from the
  IMSL Fortran library in order to perform constrained
  likelihood optimisation. The likelihood is maximised
  ensuring that the detection function remains monotonic.
  This is done by evaluating the detection function at 20
  equally-spaced distances, ensuring that the probability
  of detection is always decreasing with increasing
  distance.

  NB: this function is ONLY available under Windows.
  Covariates are not handled at the moment.

  This routine should only be called from higher level
  functions. YMMV.
}
\alias{foptim}
\alias{foptim.dll}
\value{fitted model object \item{par}{parameter estimates}
\item{fvalue}{value of the log-likelihood at the parameter estimates}
\item{message}{if the procedure has converged the message "CONVERGENCE"}
\item{aux}{list of other auxiliary information: z - covariate matrix (not
used)

zdim - dimension of z (not used)

itercept.only - (not used)

int.range - integration range (not used)

width - truncation distance

x - data frame

ftype - key function type

adj.series - adjustment series type ("NULL" for none)

adj.order - adjustment order

adj.scale - scaling for adjustments}

\item{misc.options}{contents of misc.options, above} \item{bounded}{TRUE if
the parameters hit their upper or lower bounds}}
\author{David L. Miller \email{mrds@ninepointeightone.net}}
\arguments{
  \item{initialvalues}{vector of initial values of the parameters}
  \item{bounds}{list with two elements, "upper" and "lower" containing the
upper and lower bounds for parameters.}
  \item{detfct.options}{options defining detection function; list with
elements: ftype - key function, "unif", "hn" or "hr"

adj.series - adjustment terms to use: "NULL", "poly", "cos" or "herm"

adj.order - order of adjustment

adj.scale - scaling used for the adjustment terms (scale or width)

width - truncation distance}
  \item{misc.options}{miscellaneous options; list with elements: showit - as
in other functions, shows debugging information

nonmono - should monotonicity be enforced (default no)

strict - should strict monotonicity (default no)

maxiter - maximum number of iterations for the Fortran code (default 50)

x - data frame}
}
