#' Two-part normal key function
#'
#' The two-part normal detection function of Becker and Christ (2015). Either
#' side of an estimated apex in the distance histogram has a half-normal
#' distribution, with differing scale parameters. Covariates may be included
#' but affect both sides of the function.
#'
#' Two-part normal models have 2 important parameters:
#' \itemize{
#'  \item The apex, which estimates the peak in the detection function (where
#'  g(x)=1). The log apex is reported in \code{summary} results, so taking the
#'  exponential of this value should give the peak in the plotted function (see
#'  examples).
#'  \item The parameter that controls the difference between the sides
#'  \code{.dummy_apex_side}, which is automatically added to the formula for a
#'  two-part normal model. One can add interactions with this variable as
#'  normal, but don't need to add the main effect as it will be automatically
#'  added.
#' }
#'
#' @param distance perpendicular distance vector
#' @param ddfobj meta object containing parameters, design matrices etc
#'
#' @return a vector of probabilities that the observation were detected given
#' they were at the specified distance and assuming that g(mu)=1
#'
#' @aliases two-part-normal
#' @author Earl F Becker, David L Miller
#' @references
#' Becker, E. F., & Christ, A. M. (2015). A Unimodal Model for Double Observer
#' Distance Sampling Surveys. PLOS ONE, 10(8), e0136403.
#' \doi{10.1371/journal.pone.0136403}
keyfct.tpn <- function(distance, ddfobj){

  # decide which side of the apex each observation lies on
  apex <- exp(ddfobj$shape$parameters)
  ind <- distance < apex

  # if we have just one element (e.g., when integrating) then fill this
  # out to have same # rows as distance
  xmat <- ddfobj$xmat
  if(nrow(ddfobj$scale$dm)==1){
    xmat <- xmat[rep(1, length(distance)), , drop=FALSE]
  }
  # following code convention from Earl
  xmat$.dummy_apex_side <- as.numeric(!ind)
  scale.dm <- setcov(xmat, ddfobj$scale$formula)

  # return vector
  ret <- rep(NA, length(distance))

  # find values for distances > apex
  if(sum(!ind)>0){
    right_scale <- scalevalue(ddfobj$scale$parameters,
                              scale.dm[!ind, , drop=FALSE])
    right_m <- (distance[!ind] - apex)
    right_vals <- exp( -((right_m/ (sqrt(2) * right_scale))^2) )
    ret[!ind] <- right_vals
  }

  # find values for distances <= apex
  if(sum(ind)>0){
    # get scale parameters
    left_scale <- scalevalue(ddfobj$scale$parameters,
                             scale.dm[ind, , drop=FALSE])
    # evaluate normal density
    left_m <- (distance[ind] - apex)
    left_vals <- exp( -((left_m/ (sqrt(2) * left_scale))^2) )
    # fill-in the correct values
    ret[ind] <- left_vals
  }

  return(ret)
}
