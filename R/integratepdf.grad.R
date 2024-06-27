#' Numerically integrate non-normalised pdf of observed distances over 
#' specified ranges
#' 
#' Gradient of the integral of the detection function, i.e., d beta/d theta in 
#' the documentation. This gradient of the integral is the same as the integral
#' of the gradient, thanks to Leibniz integral rule. 
#' 
#' @param par.index the index of the parameter of interest
#' @param distance vector of distances
#' @param ddfobj the ddf object
#' @param width the truncation width
#' @param int.range vector with the lower and upper bound of the integration
#' @param left the left truncation. Defaults to zero.
#' @param pdf.based evaluate the non-normalised pdf or the detection function? 
#'                  Defaults to TRUE. 
#' @param standardize not being used really, so can probably be removed.
#' 
#' @author Felix Petersma

integratepdf.grad <- function(par.index, ddfobj, int.range, width, 
                              standardize = FALSE, point = FALSE, left = 0, 
                              pdf.based = TRUE) {
  
  ## If the non-standardised detection function is required, simply 
  ## integrate the gradient of the detection function directly.
  if (!standardize) {
    out <- integrate(distpdf.grad, lower = int.range[1], upper = int.range[2],
                     par.index = par.index, ddfobj = ddfobj, width = width,
                     standardize = FALSE, point = point, left = left,
                     pdf.based = pdf.based)$value
  } 
  
  return(out)
}
