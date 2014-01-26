#' Logistic as a function of covariates
#'
#' treats logistic as a function of covariates; for a given covariate
#' combination it computes function at with those covariate values at a
#' range of distances
#'
#' @param distance vector of distance values
#' @param x covariate data
#' @param models model list
#' @param beta logistic parameters
#' @param transect \code{"line"} or \code{"point"} transect model
#'
#' @return vector of probabilities
#' @author Jeff Laake
logisticbyx <- function (distance, x, models, beta, transect){

  # Functions used: g0, setcov

  xlist <- as.list(x)
  xlist$distance <- distance
  xmat <- expand.grid(xlist)

  if(transect=="line"){
    return(g0(beta, setcov(xmat, models$g0model)$cov))
  }else{
    return(g0(beta, setcov(xmat, models$g0model)$cov)*2*distance)
  }
}
