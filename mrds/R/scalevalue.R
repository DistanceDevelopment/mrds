#' Compute scale of detection function
#'
#' Uses a log link
#'
#' @param key.scale scale parameters
#' @param z design matrix for scale covariates
#'
#' @return Vector of scale values
scalevalue <- function(key.scale, z){
  exp(as.matrix(z) %*% key.scale)
}
