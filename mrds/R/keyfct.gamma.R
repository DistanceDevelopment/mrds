#' Gamma key function
#'
#' @param distance perpendicular distance vector
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#'
#' @return vector of probabilities
keyfct.gamma <- function(distance, key.scale, key.shape){
  fr <- (1/gamma(key.shape)) * (((key.shape - 1)/exp(1))^(key.shape - 1))
  v1 <- distance/(key.scale * fr)
  return(v1^(key.shape-1)*exp(-v1)/(gamma(key.shape)*fr))
}
