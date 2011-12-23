\name{dht.se}
\alias{dht.se}
\title{Variance and confidence intervals for density and abundance estimates}
\usage{
  dht.se(model, region.table, samples, obs, options,
    numRegions, estimate.table, Nhat.by.sample)
}
\arguments{
  \item{model}{ddf model object}

  \item{region.table}{table of region values}

  \item{samples}{table of samples(replicates)}

  \item{obs}{table of observations}

  \item{options}{list of options that can be set (see
  \code{\link{dht}})}

  \item{numRegions}{number of regions}

  \item{estimate.table}{table of estimate values}

  \item{Nhat.by.sample}{estimated abundances by sample}
}
\value{
  List with 2 elements: \item{estimate.table}{completed
  table with se, cv and confidence limits} \item{vc
  }{correlation matrix of estimates}
}
\description{
  Computes standard error, cv, and log-normal confidence
  intervals for abundance and density within each region
  (if any) and for the total of all the regions.  It also
  produces the correlation matrix for regional and total
  estimates.
}
\details{
  The variance has two components: 1) variation due to
  uncertanity from estimation of the detection function and
  2) variation in abundance due to random sample selection.
  The first component is computed using a delta method
  estimate of variance (\code{\link{DeltaMethod}} (Huggins
  1989, 1991, Borchers et al. 1998) in which the first
  derivatives of the abundance estimator with respect to
  the parameters in the detection function are computed
  numerically.  The second component can be computed in one
  of three ways as set by the option \code{varflag} with
  values 0,1,2.

  A value of 0 is to use a binomial variance for the number
  of observations and it is only useful if the sampled
  region is the survey region and the objects are not
  clustered which will not occur very often. If covered
  region is less than the survey region the variance
  estimator is scaled up but it will be a poor estimator
  and the confidence interval will likely not achieve the
  nominal level.

  A value of 1 uses the variance for the encounter rate of
  (Fewster et al. (2009), estimator R2 (which has been
  shown to have better properties than the previous default
  of Buckland et al. 2001 pg 78-79)).  If
  \code{group=FALSE} the variance of the mean group size is
  also included. This variance estimator is not appropriate
  if \code{size} or a derivative of \code{size} is used in
  the any of the detection function models.

  In general if any covariates are used in the models, the
  default option 2 is preferable.  It uses a variance
  estimator based on that suggested by Innes et al. (2002)
  which used the formula for the variance ecounter rate but
  replaces the number of observations per sample with the
  estimated abundance per sample.  The difference between
  the version used here and that in Innes et al. (2002) is
  that Innes et al. use an estimator with form similar to
  that of Buckland et al. (2001), while the estimator here
  uses a form based on Fewster et al. (2009, estimator R2).

  For more on encounter rate variance estimation, see
  \link{varn}.

  Exceptions to the above occur if there is only one sample
  in a stratum. In that case it uses Poisson assumption
  (var(x)=x) and it assumes a known variance so z=1.96 is
  used for critical value. In all other cases the degrees
  of freedom for the t-distribution assumed for the
  log(abundance) or log(density) is based on the
  Satterthwaite approximation (Buckland et al. 2001 pg 90)
  for the degrees of freedom (df).  The df are weighted by
  the squared cv in combining the two sources of variation
  because of the assumed log-normal distribution because
  the components are multiplicative.  For combining df for
  the sampling variance across regions they are weighted by
  the variance because it is a sum across regions.

  A non-zero correlation between regional estimates can
  occur from using a common detection function across
  regions.  This is reflected in the correlation matrix of
  the regional and total estimates which is given in the
  value list.  It is only needed if subtotals of regional
  estimates are needed.
}
\note{
  This function is called by \code{dht} and it is not
  expected that the user will call this function directly
  but it is documented here for completeness and for anyone
  expanding the code or using this function in their own
  code
}
\author{
  Jeff Laake
}
\references{
  see \code{\link{dht}}
}
\seealso{
  \code{\link{dht}}, \code{\link{print.dht}}
}
\keyword{utility}

