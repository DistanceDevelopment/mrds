#' The gradient of the half-normal key function
#' 
#' The key function contains one parameter, the scale. 
#' Current implementation assumes that scaled dist is x/scale, not x/width
#' 
#' d key / d scale = exp(-y ^ 2 / (2 scale ^ 2)) * (y ^ 2 / scale ^ 3)
#' 
#' @param distance perpendicular distance vector
#' @param key.scale vector of scale values
#' 
#' @return vector of derivatives of the half-normal key function w.r.t. the
#' scale parameter
keyfct.grad.hn <- function(distance, key.scale){
  return(exp(-distance ^ 2 / (2 *  key.scale ^ 2)) *
             (distance ^ 2 / key.scale ^ 3))
  
  # # alternative when including the denominator \sqrt{\pi * \sigma ^2 / 2}
  # DOES NOT WORK
  # return( 
  #   -(sqrt(2) * (key.scale^2 - distance ^2) * exp(-distance ^ 2 / (2 *  key.scale ^ 2))
  #     / sqrt(pi) * key.scale^3 * abs(key.scale)))
}
