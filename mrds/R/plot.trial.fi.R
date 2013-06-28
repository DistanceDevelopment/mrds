#' Plot fit of detection functions and histograms of data from distance
#' sampling model
#' 
#' Plots the fitted detection functions for a distance sampling model and
#' histograms of the distances (for unconditional detection functions) or
#' proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data.  The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}.  The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#' 
#' In earlier versions of this package the plot routines other than
#' \code{plot.ds} (\code{plot.io} etc) used \code{plot.cond} and
#' \code{plot.uncond} to construct the conditional and unconditional detection
#' functions respectively.  This is no longer the case although the advanced
#' user can still call \code{plot.cond} and \code{plot.uncond} if required.  It
#' is not intended for the user to call any of \code{plot.ds},
#' \code{plot.trial.fi}, \code{plot.trial},\code{plot.rem.fi}, \code{plot.rem},
#' \code{plot.io.fi} or \code{plot.io} but the arguments are documented here.
#' Instead the generic \code{plot} command should be used and it will call the
#' appropriate function based on the type of \code{ddf} object.
#' 
#' For plot routine \code{plot.ds} the \code{which} command allows the user to
#' select which plots are returned from the following options: \tabular{ll}{
#' \code{which} \tab \code{plot}\cr \code{1} \tab data summary plot - a
#' histogram of the observed distances \cr \code{2} \tab a scaled histogram of
#' detections with a line giving the detection function averaged over the
#' estimated population levels of the covariate values, and one dot for each
#' observation at its estimated detection probability.\cr }
#' 
#' For plot routines \code{plot.trial.fi} and \code{plot.trial} the following
#' plots are available: \tabular{ll}{ \code{which} \tab \code{plot} \cr
#' \code{1} \tab data summary plot - a histogram of the observed distances for
#' observer 1 \cr \code{2} \tab data summary plot - a histogram of the observed
#' distances for observer 2 \cr \code{3} \tab Observer 1 detection function - a
#' scaled histogram of detections with fitted DS model scaled from the MR
#' estimated g(0).  The line shows the population average detection function
#' and the points display estimated detection probability \cr \code{4} \tab
#' Conditional MR detection function - observer 1 given obs 2, giving the
#' proportion of duplicates with fitted MR model averaged over population
#' covariate vales and dots for each estimated detection probability. \cr } For
#' plot routines \code{plot.rem.fi} and \code{plot.rem} the following plots are
#' available: \tabular{ll}{ \code{which} \tab \code{plot} \cr \code{1} \tab
#' data summary plot - a histogram of the observed distances for observer 1 \cr
#' \code{2} \tab data summary plot - a histogram of the observed distances for
#' observer 2 \cr \code{3} \tab Observer 1 detection function - a scaled
#' histogram of detections with fitted DS model scaled from the MR estimated
#' g(0).  The line shows the population average detection function and the
#' points display estimated detection probability \cr \code{4} \tab Conditional
#' MR detection function - observer 1 given obs 2, giving the proportion of
#' duplicates with fitted MR model averaged over population covariate vales and
#' dots for each estimated detection probability. \cr }
#' 
#' For plot routines \code{plot.io.fi} and \code{plot.io} the following plots
#' are available: \tabular{ll}{ \code{which} \tab \code{plot} \cr \code{1} \tab
#' data summary plot - a histogram of the observed distances for observer 1 \cr
#' \code{2} \tab data summary plot - a histogram of the observed distances for
#' observer 2 \cr \code{3} \tab Observer 1 detection function - a scaled
#' histogram of detections with fitted DS model scaled from the MR estimated
#' g(0).  The line shows the population average detection function and the
#' points display estimated detection probability \cr \code{4} \tab Observer 2
#' detection function - as for \code{plot} 3 but using the detections from
#' Observer 2\cr \code{5} \tab Duplicates detection function - as for
#' \code{plot} 3 but using the duplicate detections\cr \code{6} \tab Pooled
#' detection function - as for \code{plot} 3 but using the pooled detections\cr
#' \code{7} \tab Conditional MR detection function - observer 1 given obs 2,
#' giving the proportion of duplicates with fitted MR model averaged over
#' population covariate vales and dots for each estimated detection
#' probability. \cr \code{8} \tab Conditional MR detection function - observer
#' 2 given obs 1, giving the proportion of duplicates with fitted MR model
#' averaged over population covariate vales and dots for each estimated
#' detection probability. \cr }
#' 
#' @aliases plot.trial.fi
#' @method plot trial.fi
#' @S3method plot trial.fi
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
#' @param xlab label for x-axis
#' @param ylab label for y-axis 
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (plot, hist, lines, points, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot


"plot.trial.fi" <-
		function(x, which=1:2, breaks=NULL, nc=NULL,  maintitle="", showlines=TRUE, showpoints=TRUE, 
				ylim=c(0,1),angle=-45,density=20,col="black",jitter=NULL,divisions=25,xlab="Distance",ylab="Detection probability",...)
#
# plot.trial.fi
#
# Provides plots of fitted functions for a trial.fi object
#
#  Functions used: process.data, predict(predict.trial.fi),plot.uncond, plot.cond, calcp.mrds
#
{
#
#  Retrieve values from model object
#
	model <- x
	data<- process.data(model$data,model$meta.data)$xmat
	xmat<- data[data$observer==1&data$detected==1,]
	gxvalues <- predict(model,newdata=xmat,type="response",integrate=FALSE)$fitted
	width <- model$meta.data$width 
	left<- model$meta.data$left
#
#  If number of classes for histogram intervals was not set compute a reasonable default
#
	if(is.null(nc))
		nc<-round( sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
								length(xmat$distance[xmat$observer==1&xmat$timesdetected==2]) )),0)
	hascov <- TRUE
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
		if(.Platform$GUI=="Rgui")dev.new()
	   plot_uncond(model,1,xmat,gxvalues,nc,finebr=(width/divisions)*(0:divisions),breaks,showpoints,showlines,
			 maintitle,ylim,point=model$meta.data$point,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
   }
#
#  Plot conditional detection functions
#
    if(is.element(2,which))
	{
		if(.Platform$GUI=="Rgui")dev.new()
	   xmat<- data[data$observer==2&data$detected==1,]
	   gxvalues <- predict(model,newdata=xmat,type="response",integrate=FALSE)$fitted  
	   est<-calcp.mrds(model$mr$formula,model$mr$family$link,model$mr$coefficients,model$data,vname="distance",
			        lower=left,upper=width,divisions=divisions,type=model$meta.data$point,objname="object",obsname="observer")
	   plot_cond(1,data,gxvalues,list(x=est$x,p=est$p1),nc,breaks,showpoints,showlines,
			maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,ylab=ylab,...)
    }
	invisible(NULL)
}


