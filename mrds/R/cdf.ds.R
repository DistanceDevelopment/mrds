


#' Cumulative distribution function (cdf) for fitted distance sampling
#' detection function
#' Computes cdf values of observed distances from fitted distribution.  For a
#' set of observed x it returns the integral of f(x) for the range= (inner, x),
#' where inner is the innermost distance which is observable (either 0 or left
#' if left truncated).  In terms of g(x) this is the integral of g(x) over
#' range divided by the integral of g(x) over the entire range of the data
#' (inner, W).
#' 
#' 
#' @param model fitted distance sampling model
#' @param newdata new data values if computed for values other than the
#'   original observations
#' @return vector of cdf values for each observation
#' @note This is an internal function that is not intended to be invoked
#'   directly.  It is called by \code{\link{qqplot.ddf}} to compute values for
#'   K-S and CvM tests and the Q-Q plot.
#' @author Jeff Laake
#' @seealso \code{\link{qqplot.ddf}}
#' @keywords utility
cdf.ds <-
function(model,newdata=NULL)
{
#
# cdf.ds - Computes cdf values of observed distances from fitted distribution.  For a set of observed x
# it returns the integral of f(x) for the range= (inner, x), where inner is the innermost distance which is
# observable (either 0 or left if left truncated).  In terms of g(x) this is the integral of g(x) over range divided
# by the integral of g(x) over the entire range of the data (inner, W).  the
# is somewhat more complicated to allow for variable integration ranges (int.range).
#
# arguments:
#
# model    - object from fit.ds
# newdata  - new data for computation if any
#
#
# return value
#  cdf for each detection
# 
# Uses: integratedetfct, integratedetfct.logistic
#   
# 
#      Extract values from model
#
  ltmodel <- model$ds
  fpar <- model$par
  width <- ltmodel$aux$width
  ddfobj <-ltmodel$aux$ddfobj 
  x <- ddfobj$xmat  
  z <- ddfobj$scale$dm
  zdim <- dim(z)
  ftype <- ddfobj$type
  intercept.only <- ddfobj$intercept.only
  cgftab <- ddfobj$cgftab
#
# Set up integration range
#
  if(is.null(ltmodel$aux$int.range)){
    int.range=cbind(rep(0,dim(x)[1]+1),c(width,x$distance))
  }else{
    int.range=ltmodel$aux$int.range
    if(is.matrix(int.range))
      int.range[,2]=c(width,x$distance)
    else
      int.range=cbind(rep(int.range[1],dim(x)[1]+1),c(width,x$distance))
  }
#
#      Do integration of g(x) from inner (0 or left) to x
#
#  dlm 11/07/2005 Commented out logistic stuff
#
#       if(ftype=="logistic")
#          int1=integratedetfct.logistic (x,ltmodel$model$scalemodel,width,int.range,theta1,ltmodel$aux$integral.numeric,z)
#       else

  if(any(is.null(ddfobj$adjustment$order))){
	  int1 <- integratedetfct(ddfobj=ddfobj,select=rep(TRUE,nrow(x)),width=width,int.range=int.range,
			  doeachint=TRUE,standardize=ltmodel$aux$misc.options$standardize,point=FALSE)
	  
#    int1=integratedetfct(cgftab,ftype,width, int.range, z, NULL, fpar,fpar,intercept.only, FALSE,
#                         adj.series=adj.series,adj.order=adj.order,adj.scale=adj.scale)
#
#   Divide by integral of g(x) over entire integration range (e.g., 0 to W); thus providing
#   integral of f(x) from inner to x.
#
    int2=predict(model,integrate=TRUE,esw=FALSE,compute=TRUE)$fitted
		
# dlm 29-Aug-05	Trying a different method...

#  int2=integratedetfct(cgftab,ftype,width, int.range=c(0,width), z, NULL,fpar,fpar,
#                              intercept.only,FALSE,adj.series=adj.series,adj.order=adj.order,
#                              adj.scale=adj.scale,standardize=FALSE)

  }else{

    left <- int.range[2:dim(int.range)[1],1]
    right <- int.range[2:dim(int.range)[1],2]

    int1 <- integrate(detfct,select=rep(TRUE,nrow(x)),lower=left,upper=right,ddfobj=ddfobj,
                      width=width)$value

	int2 <- integrate(detfct,select=rep(TRUE,nrow(x)),lower=0,upper=width,ddfobj=ddfobj,
					  width=width)$value

  }

  fitted=int1/int2

  return(list(fitted=fitted))
}
