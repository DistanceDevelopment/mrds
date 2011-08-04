#' Integral of standardized (scale=1) detection function
#' 
#' Computes the integral of a standardized (scale=1) detection function for a
#' specified range
#' 
#' If \code{xmin} is specified it computes the integral from \code{xmin} to
#' each value in \code{x}; otherwise, it computes the integral from 0 to
#' \code{x}.
#' 
#' @param x vector of lower or upper values for integral bound (depending on
#'   value of xmin)
#' @param xmin lower bound if specified
#' @param ddfobj distance detection function specification
#' @param index specific data row index
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count(TRUE) or line
#'   transect(FALSE)
#' @return vector of integral values of detection function with scale=1
#' @note This is an internal function that is not intended to be invoked
#'   directly.
#' @author Jeff Laake
#' @seealso \code{\link{tablecgf}}
#' @keywords utility
gstdint <-
function (x, xmin = NULL, ddfobj, index=NULL,select=NULL,width,standardize=TRUE, point=FALSE)
#
#  gstdint - computes integrals of standardized detection function (scale=1)
#
#  Arguments:
#  
# x             - vector of lower or upper values for integral bound (depending on value of xmin) 
# xmin          - lower bound if specified
# ddfobj        -distance detection function specification
# index         -specific data row index
# select        -logical vector for selection of data values
# width         -truncation width
# standardize   -logical used to decide whether to divide through by the function evaluated at 0 }
# point         -logical to determine if point count(TRUE) or line transect(FALSE)
#
#  Value: vector of integral values of det fct with scale=1
#
#  17-Aug-05; This function was changed to use unequal set of grid points and a reduced
#  integration tolerance to speed up execution.
#
#  21-Aug-05; dlm - changed code so that it works with new structure
#
#  19-Jan-06; jll modified for new detfct so it would std integral with scale=1; stdint=TRUE
#  2 Sept-10; now works for points
{
  	if(!point)
	{
  		if(is.null(xmin)){
    		return(integrate(fx, lower = x[1], upper = x[2], width=width,
           		ddfobj=ddfobj, select=select, index=index, rel.tol = 1e-7,standardize=standardize,
           		stdint=TRUE)$value)
  		}else{ 
    		return(integrate(fx, lower = xmin, upper = x, width=width,
		   		ddfobj=ddfobj, select=select, index=index, rel.tol = 1e-7,standardize=standardize,
           		stdint=TRUE)$value)	
		}
	}
	else
	{
		if(is.null(xmin)){
			return(integrate(fr, lower = x[1], upper = x[2], width=width,
							ddfobj=ddfobj, select=select, index=index, rel.tol = 1e-7,standardize=standardize,
							stdint=TRUE)$value)
		}else{ 
			return(integrate(fr, lower = xmin, upper = x, width=width,
							ddfobj=ddfobj, select=select, index=index, rel.tol = 1e-7,standardize=standardize,
							stdint=TRUE)$value)				
		}
  	}
}

