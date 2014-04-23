\name{keyfct.hz}
\alias{keyfct.hz}
\title{Hazard rate key function}
\usage{
keyfct.hz(distance, key.scale, key.shape)
}
\arguments{
  \item{distance}{perpendicular distance vector}

  \item{key.scale}{vector of scale values}

  \item{key.shape}{vector of shape values}
}
\value{
vector of probabilities that the observations were detected
given they were at the specified distance and assuming that
g(0)=1 (ie a standard line transect detection function).
}
\description{
Function has the form 1 - exp (- (x/scale) ^ (-shape)))
}

