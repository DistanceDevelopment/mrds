\name{ddf.io.fi}
\alias{ddf.io.fi}
\title{Mark-Recapture Analysis of Independent Observer Configuration with Full...}
\usage{\method{ddf}{io.fi} (model, data, meta.data=list(),
    control=list(), call="", method)}

\description{
  Mark-Recapture Analysis of Independent Observer
  Configuration with Full Independence The mark-recapture
  data derived from an independent observer distance
  sampling survey can be used to derive conditional
  detection functions (p_j(y)) for both observers (j=1,2).
  They are conditional detection functions because
  detection probability for observer j is based on seeing
  or not seeing observations made by observer 3-j. Thus,
  p_1(y) is estimated by p_1|2(y).  If detections by the
  observers are independent (full independence) then
  p_1(y)=p_1|2(y),p_2(y)=p_2|1(y) and for the union, full
  independence means that p(y)=p_1(y) + p_2(y) -
  p_1(y)*p_2(y) for each distance y.  In fitting the
  detection functions the likelihood given by eq 6.8 and
  6.16 in Laake and Borchers (2004) is used. That analysis
  does not require the usual distance sampling assumption
  that perpendicular distances are uniformly distributed
  based on line placement that is random relative to animal
  distribution.  However, that assumption is used in
  computing predicted detection probability which is
  averaged based on a uniform distribution (see eq 6.11 of
  Laake and Borchers 2004).
}

\details{
  For a complete description of each of the calling
  arguments, see \code{\link{ddf}}.  The argument
  \code{model} in this function is the same as
  \code{mrmodel} in \code{ddf}.  The argument
  \code{dataname} is the name of the dataframe specified by
  the argument \code{data} in \code{ddf}. The arguments
  \code{control},\code{meta.data},and \code{method} are
  defined the same as in \code{ddf}.
}
\value{result: an io.fi model object}
\author{Jeff Laake}
\seealso{\code{\link{ddf.io}},\code{\link{summary.io.fi}},\code{\link{coef.io.fi}},\code{\link{plot.io.fi}},
\code{\link{gof.io.fi}},\code{\link{io.glm}}}
\references{Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
Thomas. Oxford University Press.}
\keyword{Statistical}
\keyword{Models}
\arguments{
  \item{model}{mark-recapture model specification}
  \item{data}{analysis dataframe}
  \item{meta.data}{list containing settings controlling data structure}
  \item{control}{list containing settings controlling model fitting}
  \item{call}{original function call used to call \code{ddf}}
  \item{method}{analysis method; only needed if this function called from
\code{ddf.io}}
}
