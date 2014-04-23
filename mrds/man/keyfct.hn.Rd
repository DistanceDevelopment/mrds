\name{keyfct.hn}
\alias{keyfct.hn}
\title{Half-normal key function
Function has the form exp (- (x^2/(2*sigma^2))
The 2 is included in the denominator to preserve the
relationship with the standard deviation in a normal distribution.}
\usage{
keyfct.hn(distance, key.scale)
}
\arguments{
  \item{distance}{perpendicular distance vector}

  \item{key.scale}{vector of scale values}
}
\value{
a vector of probabilities that the observation were
detected given they were at the specified distance and
assuming that g`(0)=1 (ie a standard line transect
detection function)
}
\description{
Half-normal key function Function has the form exp (-
(x^2/(2*sigma^2)) The 2 is included in the denominator to
preserve the relationship with the standard deviation in a
normal distribution.
}

