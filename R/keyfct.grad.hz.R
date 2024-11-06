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
#' When distance = 0, the gradients are also zero. However, the equation below
#' will result in NaN and (-)Inf due to operations such as log(0) or division by
#' zero. We correct for this in line 33.
#' 
#' @param distance perpendicular distance vector
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#' @param shape is the gradient parameter the shape parameter? Defaults to FALSE
#' 
#' @return matrix of derivatives of the hazard rate key function w.r.t. the
#' scale parameter and the shape parameter.
keyfct.grad.hz <- function(distance, key.scale, key.shape, shape = FALSE){
if (!shape) {
  out <- (key.shape * exp(-(1/ (distance / key.scale) ^ key.shape) ) ) / 
           ( ( (distance / key.scale) ^ key.shape ) * key.scale) 
} else {
  out <- -1 * ((log(distance / key.scale) * 
                  exp(-(1/ (distance / key.scale) ^ key.shape))) /
                 (distance / key.scale) ^ key.shape)
}
  ## Set out to zero where distance is zero
  out[distance == 0] <- 0
  
  ## Return the result
  return(out)
}
# keyfct.grad.hz <- function(distance, key.scale, key.shape){
#   out.scale <- (key.shape * exp(-(1/ (distance / key.scale) ^ key.shape)) / 
#                   ((distance / key.scale) ^ key.shape ) * key.scale)
#   out.shape <- -1 * ((log(distance / key.scale) * 
#                         exp(-(1/ (distance / key.scale) ^ key.shape))) /
#                        (distance / key.scale) ^ key.shape)
#   return(matrix(c(out.scale, out.shape), ncol = 2, byrow = FALSE))
# }
