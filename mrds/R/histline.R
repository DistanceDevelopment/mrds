#' Plot histogram line
#' 
#' Takes bar heights (height) and cutpoints (breaks), and constructs a
#' line-only histogram from them using the function plot() (if lineonly==FALSE)
#' or lines() (if lineonly==TRUE).
#' 
#' 
#' @param height heights of histogram bars
#' @param breaks cutpoints for x
#' @param lineonly if TRUE, drawn with plot; otherwise with lines to allow
#'   addition of current plot
#' @param outline if TRUE, only outline of histogram is plotted
#' @param fill If fill==TRUE, uses polygon() to fill bars
#' @param ylim limits for y axis
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param det.plot if TRUE, plot is of detection so yaxis limited to unit
#'   interval
#' @param \dots Additional unspecified arguments for plot (fill==TRUE
#' @return None
#' @author ???
histline<-function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),xlab="x",ylab="y",det.plot=FALSE,...)
#-------------------------------------------------------------------------------------
# Takes bar heights (height) and cutpoints (breaks), and constructs a line-only
# histogram from them using the function plot() (if lineonly==FALSE) or lines()
# (if lineonly==TRUE).
# If fill==TRUE, uses polygon() to fill bars
# If fill==TRUE, valid arguments to plot() or lines() are passed via argument(s) "..."
# If outline==TRUE, only outline of histogram is plotted
# If fill!=TRUE, valid arguments to polygon() are passed via argument(s) "..."
#
#-------------------------------------------------------------------------------------
{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      if(det.plot){plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab,yaxp=c(0,1,5))}
      else{plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)}
      polygon(x,y,...)
    }
  }
}
