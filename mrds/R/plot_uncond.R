#' Plot unconditional detection function from distance sampling model
#' 
#' Plots unconditional detection function for observer=obs observations 
#' overlays histrogram, average detection function and values for individual observations
#' data.  
#'  
#' @aliases plot_uncond
#' @param model   fitted model from \code{ddf}
#' @param obs value of observer for plot
#' @param xmat processed data
#' @param gxvalues detection function values for each observation
#' @param nc number of equal-width bins for histogram
#' @param finebr fine break values over which line is averaged
#' @param breaks user define breakpoints
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE plots average predicted value line
#' @param maintitle main title line for each plot
#' @param ylim range of y axis; defaults to (0,1)
#' @param point point count data if TRUE
#' @param return.lines if TRUE, returns values for line
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour 
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter. 
#' @param xlab label for x-axis
#' @param ylab label for y-axis 
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (plot, hist, lines, points, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot
#' @examples
#' \donttest{
#' data(book.tee.data)
#' region<<-book.tee.data$book.tee.region
#' egdata<<-book.tee.data$book.tee.dataframe
#' samples<<-book.tee.data$book.tee.samples
#' obs<<-book.tee.data$book.tee.obs
#' xx=ddf(dsmodel = ~mcds(key = "hn", formula = ~sex), data = egdata[egdata$observer==1, ], 
#'  method = "ds", meta.data = list(width = 4))
#' plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1)
#' }
"plot_uncond" <-
function(model,obs,xmat,gxvalues,nc,finebr,breaks,showpoints,showlines,maintitle,ylim,point=FALSE,return.lines=FALSE,angle=-45,density=20,col="black",jitter=NULL,
		  xlab="Distance",ylab="Detection probability",...)
#
#  Value: if(return.lines) returns dataframe from average.line, otherwise NULL
#
#  Functions Used: average.line
#
{
  if(obs<=2)
  { 
     n <- length(xmat$distance[xmat$observer == obs&xmat$detected==1])
     selmat <- xmat[xmat$observer== obs,]
     det.detected <- xmat$detected[xmat$observer== obs]==1
  } else
  if(obs==3)
  {   
    n <- length(xmat$distance[xmat$observer==1])
    selmat <- xmat[xmat$observer==1,]
    det.detected<- selmat$observer==1
  } else
  {
     n <- length(xmat$distance[xmat$observer==1 & xmat$timesdetected==2])
     selmat <- xmat[xmat$observer==1,]
     det.detected <- selmat$timesdetected==2
  }
  hist.obj<-hist(selmat$distance[det.detected], breaks=breaks, plot = FALSE)
  if(!point)
	  expected.counts<-((breaks[2:(nc+1)]-breaks[1:nc])*(model$Nhat/breaks[nc+1] ))
  else
	  expected.counts<--apply(matrix(c(breaks[2:(nc+1)]^2,breaks[1:nc]^2),ncol=2,nrow=nc),1,diff)*(model$Nhat/breaks[nc+1]^2 )
  hist.obj$density<-hist.obj$counts/(expected.counts)
  hist.obj$intensities<-hist.obj$density
  freq<-hist.obj$density
  hist.obj$equidist<-FALSE
  line=average.line(finebr,obs,model)
  linevalues <- line$values
  xgrid <- line$xgrid
  maxl <- max(freq, gxvalues)
  ylim<-c(0,max(ylim,hist.obj$density)) # *** DLB change ***
  histline(hist.obj$density,breaks=breaks,lineonly=FALSE,xlab=xlab,ylab=ylab,
		     ylim=ylim,fill=TRUE,angle=angle,density=density,col=col,det.plot=TRUE,...)
  if(showlines)lines(xgrid,linevalues,...)
  if(showpoints)
  {
	  ifelse(is.null(jitter),jitter.p<-1,jitter.p<-rnorm(length(gxvalues),1,jitter))
	  points(selmat$distance,gxvalues*jitter.p,...)
  }
  if(maintitle!="")
     if(obs<=2)
        title(paste(maintitle, "\nObserver = ",obs, " detections"),...)
     else
        if(obs==3)
           title(paste(maintitle, "\nPooled detections"),...)
        else
           title(paste(maintitle, "\nDuplicate detections"),...)
 if(return.lines) invisible(line) else invisible(NULL) # *** dlb change ***
}
