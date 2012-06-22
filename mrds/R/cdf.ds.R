#' Cumulative distribution function (cdf) for fitted distance sampling
#' detection function
#' 
#' Computes cdf values of observed distances from fitted distribution.  For a
#' set of observed x it returns the integral of f(x) for the range= (inner, x),
#' where inner is the innermost distance which is observable (either 0 or left
#' if left truncated).  In terms of g(x) this is the integral of g(x) over
#' range divided by the integral of g(x) over the entire range of the data
#' (inner, W).
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
  point <- model$meta.data$point
#
# Set up integration ranges
#
  if(is.null(ltmodel$aux$int.range)){
    int.range=as.matrix(cbind(rep(0,nrow(x)+1),c(width,x$distance)))
  }else{
    int.range=ltmodel$aux$int.range
	if(is.vector(int.range))int.range=matrix(int.range,nrow=1)
    if(nrow(int.range)>1)
      int.range[,2]=c(width,x$distance)
    else
      int.range=cbind(rep(int.range[1],nrow(x)+1),c(width,x$distance))
  }
  int.range <- int.range[-1,]

#
# If there are no adjustments
#
#  if(any(is.null(ddfobj$adjustment$order))){
#
#      Do integration of g(x) from inner (0 or left) to x
#
	int1<-integratepdf(ddfobj=ddfobj,select=rep(TRUE,nrow(x)),width=width,
                          int.range=int.range, doeachint=TRUE,
                          standardize=TRUE,point=point)
	  
#   integral of g(x) over entire integration range 
    int2=predict(model,integrate=TRUE,esw=FALSE,compute=TRUE)$fitted
		
#
# If there are adjustments...
#
#  }else{

#    left <- int.range[2:dim(int.range)[1],1]
#    right <- int.range[2:dim(int.range)[1],2]

# dlm Oct-11 this doesn't work with Jeff's new detfct...
#    int1<-integrate(detfct,select=rep(TRUE,nrow(x)),lower=left,upper=right,
#                    ddfobj=ddfobj,width=width)$value
#
#    int2<-integrate(detfct,select=rep(TRUE,nrow(x)),lower=0,upper=width,
#                    ddfobj=ddfobj,width=width)$value

    # have to build this bit-by-bit

# Jeff doesn't like this
#    int1<-int2<-double(nrow(x))
#    sel<-rep(FALSE,nrow(x))
#
#    for(i in 1:length(x$distance)){
#
#      this.select<-sel
#      this.select[i]<-TRUE
#
#      int1[i]<-integrate(detfct,select=this.select,
#                         lower=left[i],upper=right[i],
#                         ddfobj=ddfobj,width=width)$value
#      int2[i]<-integrate(detfct,select=this.select,
#                         lower=0,upper=width,
#                         ddfobj=ddfobj,width=width)$value
#    }

   # this might be better?
#   n<-length(x$distance)

   # create a matrix to apply()
   # left limit, right limit, index (for select), number of entries
#   lrdata<-cbind(left,right,1:n,rep(n,n))
   
   # integration function
#  intf<-function(lrdata,width,ddfobj){
#      sel<-rep(FALSE,lrdata[4])
#      sel[lrdata[3]]<-TRUE
#      integrate(detfct,select=sel,
#                lower=lrdata[1],upper=lrdata[2],
#                ddfobj=ddfobj,width=width)$value
#   }
   
   # apply number 1
#   int1<-apply(lrdata,1,intf,width=width,ddfobj=ddfobj)
   # modify the data so it's just integration over (0,width)
#   lrdata<-cbind(rep(0,n),rep(width,n),1:n,rep(n,n))
   # do the integration 
#   int2<-apply(lrdata,1,intf,width=width,ddfobj=ddfobj)

#  }

#
#   Divide by integral of g(x) over entire integration range (e.g., 0 to W); 
#   thus providing integral of f(x) from inner to x.
#
  fitted=int1/int2

  return(list(fitted=fitted))
}
