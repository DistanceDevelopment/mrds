#' Integral of pdf of distances
#'
#' Computes the integral of distpdf with scale=1 (stdint=TRUE) or specified
#' scale (stdint=FALSE)
#'
#' @param x lower,upper value for integration
#' @param ddfobj distance detection function specification
#' @param index specific data row index
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count(TRUE) or line
#'   transect(FALSE)
#' @param stdint if TRUE, scale=1 otherwise specified scale used
#' @return vector of integral values of detection function
#' @note This is an internal function that is not intended to be invoked
#'   directly.
#' @author Jeff Laake
#' @seealso \code{\link{tablecgf}}
#' @keywords utility
gstdint <- function (x, ddfobj, index=NULL,select=NULL,width,
                      standardize=TRUE, point=FALSE, stdint=TRUE){

   return(integrate(distpdf, lower = x[1], upper = x[2], width=width,
                    ddfobj=ddfobj, select=select, index=index,
                    rel.tol = 1e-7,standardize=standardize,
                    stdint=stdint,point=point)$value)
}
