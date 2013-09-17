#' Plot fit of detection functions and histograms of data from distance
#' sampling independent observer model
#' 
#' Plots the fitted detection functions for a distance sampling model and
#' histograms of the distances (for unconditional detection functions) or
#' proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data.  
#' 
#' The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}.  The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#' 
#' It is not intended for the user to call any of \code{plot.ds},
#' \code{plot.trial.fi}, \code{plot.trial},\code{plot.rem.fi}, \code{plot.rem},
#' \code{plot.io.fi} or \code{plot.io} but the arguments are documented here.
#' Instead the generic \code{plot} command should be used and it will call the
#' appropriate function based on the type of \code{ddf} object.
#' 
#' The \code{which} command allows the user to
#' select which plots are returned. See which argument definition. 
#' 
#' @aliases plot.io.fi
#' @S3method plot io.fi
#' @method plot io.fi
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced. 1: uncond det fct, 2:cond det fct
#' @param breaks user define breakpoints
#' @param nc number of equal-width bins for histogram
#' @param maintitle main title line for each plot
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE a line representing the average
#'   detection probability is plotted
#' @param ylim range of y axis; defaults to (0,1)
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour 
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.  
#' @param divisions number of divisions for averaging line values; default = 25 
#' @param new if TRUE, opens new device for each plot; set new=FALSE if you use par(mfrow=..) or layout
#' @param xlab label for x-axis
#' @param ylab label for y-axis 
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (plot, hist, lines, points, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot

"plot.io.fi" <-
		function(x, which=1:6, breaks=NULL, nc=NULL,  maintitle="", showlines=TRUE, showpoints=TRUE, 
				ylim=c(0,1),angle=-45,density=20,col="black",jitter=NULL,divisions=25,new=TRUE,xlab="Distance",ylab="Detection probability",subtitle=TRUE,...)
#
# plot.io.fi
#
# Provides plots of fitted functions for a io.fi object
#
#  Functions used: process.data, predict(predict.io.fi), plot.uncond, plot.cond, calcp.mrds
#
{
#
#  Retrieve values from model object
#
	model <- x
	xmat<- model$mr$data
	xmat$offsetvalue=0
	cond.det <- predict(model,newdata=xmat,integrate=FALSE)  
	fitted <-cond.det$fitted
	p1<- cond.det$p1
	p2<- cond.det$p2    
	width <- model$meta.data$width 
	left <- model$meta.data$left
#
#   If number of classes for histogram intervals was not set compute a reasonable default
#
	if(is.null(nc))
		nc<-round(sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
								length(xmat$distance[xmat$observer==2&xmat$detected==1]),length(xmat$distance[xmat$observer==1&xmat$timesdetected==2]) )),0)
#
#  Set up default break points unless specified
#
	if(model$meta.data$binned)
	{
		breaks<-model$meta.data$breaks
		nc<-length(breaks)-1
	} else
	if(is.null(breaks))
		breaks <- left + ((width-left)/nc)*(0:nc)
	else
		nc=length(breaks)-1
#
#  Plot primary unconditional detection function
#
	if(is.element(1,which))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_uncond(model,1,xmat,gxvalues=p1,nc,finebr=(width/divisions)*(0:divisions),breaks,showpoints,showlines,
				maintitle,ylim,point=model$meta.data$point,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
	}
#
#  Plot secondary unconditional detection function
#
	if(is.element(2,which))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_uncond(model,2,xmat,gxvalues=p2,nc,finebr=(width/divisions)*(0:divisions),breaks,showpoints,showlines,
				maintitle,ylim,point=model$meta.data$point,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
	}
#
#  Plot pooled unconditional detection function
#
	if(is.element(3,which))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_uncond(model,3,xmat,gxvalues=p1+p2-p1*p2,nc,finebr=(width/divisions)*(0:divisions),breaks,showpoints,showlines,
				maintitle,ylim,point=model$meta.data$point,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
	}
#
#  Plot duplicate unconditional detection function
#
	if(is.element(4,which))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_uncond(model,4,xmat,gxvalues=p1*p2,nc,finebr=(width/divisions)*(0:divisions),breaks,showpoints,showlines,
				maintitle,ylim,point=model$meta.data$point,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
	}
#
#  Plot conditional detection functions
#
	data=model$mr$data
	data$offsetvalue=0
    if(is.element(5,which))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		gxvalues <-p1[xmat$detected[xmat$observer==2]==1] 
		plot_cond(1,data,gxvalues,model,nc,breaks,finebr=(width/divisions)*(0:divisions),showpoints,showlines,
				maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,subtitle=subtitle,...)
	}
#
#  Plot secondary conditional detection function
#
	if(is.element(6,which))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		gxvalues <-p2[xmat$detected[xmat$observer==1]==1] 
		plot_cond(2,data,gxvalues,model,nc,breaks,finebr=(width/divisions)*(0:divisions),showpoints,showlines,
				maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,subtitle=subtitle,...)
	}
	invisible(NULL)
}
