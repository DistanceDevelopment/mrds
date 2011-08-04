#' Spline approximation to scale dependent integration of detection function
#' 
#' Computes spline approximation to cumulative detection function for
#' interpolation of scale dependent integrations of the detection function.
#' 
#' This is an internal function used to speed up integration of the detection
#' function with scale covariates. The detection function is integrated at a
#' series of points from 0 to W and then a spline is fitted to the computed
#' values which are cumulative (integral from 0 to x < integral from 0 to x+dx
#' fr dx>0).  The spline is then used to predict values of the integral which
#' depend on the scale which can depend on observation specific covariates.
#' 
#' @param ddfobj distance detection function specification
#' @param width truncation width
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count(TRUE) or line
#'   transect(FALSE)
#' @return a smooth.spline result object
#' @note This is an internal function that is not intended to be invoked
#'   directly.
#' @author Jeff Laake
#' @seealso \code{\link{gstdint}}
#' @keywords utility
tablecgf <-
function(ddfobj, width ,standardize=TRUE,point=FALSE)
{
#
# Tablecgf 
#
# Computes spline approximation to cumulative detection function for interpolation of scale dependent integrations 
# of the detection function.
#
# Arguments:
# ddfobj         - distance sampling object
# width          - scalar, vector or matrix of integration bounds
# standardize    - logical used to decide whether to divide through by the function evaluated at 0 
# point          - logical to determine if point count(TRUE) or line transect(FALSE)
#	
#
# Functions:
# gstdint        -computes integral of standardized detection function over interval x,x+gridint
# cumsum         -splus function that returns the cumulative sums of a vector
# smooth.spline  -splus function that creates spline approximation for interpolation
# apply          -splus function which applies a named function to the elements of a vector
#
# This function was completely changed in Aug 05 to use a fixed set of grid points
# which seems to be more reliable and possibly a little quicker
#
#Create grid over interval with exponentially spaced standard set of 100 grid points
    xx=exp(0.05*(1:100))-1
    xx=cbind(c(0,xx[1:99]),xx)
#Create cumulative sums of values of g(x) integrals from grid (xx)
	y <- cumsum(apply(xx, 1, gstdint, ddfobj=ddfobj,index=1,width=width,standardize=standardize,point=point))
#Return smoothed spline of cumulative integral values
    spp <- smooth.spline(c(0,xx[,2]), c(0, y))
    return(spp)
}


