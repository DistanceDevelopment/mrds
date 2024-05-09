#' The gradient of the hazard-rate key function
#' 
#' The key function contains two parameters, the scale and the shape, and 
#' so the gradient is two-dimensional. Current implementation assumes that 
#' scaled dist is x/scale, not x/width
#' 
#' d key / d scale = (shape * exp(-(1/ (x/scale) ^ shape)) / 
#'                      ((x/scale) ^ shape ) * scale)
#' d key / d shape = - ((log(x / scale) * exp(-(1/ (x/scale) ^ shape))) / 
#'                      (x/scale) ^ shape)
#' 
#' @param distance perpendicular distance vector
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#' 
#' @return matrix of derivatives of the hazard rate key function w.r.t. the
#' scale parameter and the shape parameter.
keyfct.grad.hz <- function(distance, key.scale, key.shape){
  out.scale <- (key.shape * exp(-(1/ (distance / key.scale) ^ key.shape)) / 
                  ((distance / key.scale) ^ key.shape ) * key.scale)
  out.shape <- -1 * ((log(distance / key.scale) * 
                        exp(-(1/ (distance / key.scale) ^ key.shape))) /
                       (distance / key.scale) ^ key.shape)
  return(matrix(c(out.scale, out.shape), ncol = 2, byrow = FALSE))
}
