#' Plot conditional detection function from distance sampling model
#' 
#' Plot proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data.  
#' 
#' 
#' @aliases plot_cond
#' @param obs obsever code
#' @param xmat processed data
#' @param gxvalues detection function values for each observation
#' @param est line values constructed with calcp.mrds
#' @param nc number of equal-width bins for histogram
#' @param breaks user define breakpoints
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE plots average predicted value line
#' @param maintitle main title line for each plot
#' @param ylim range of y axis; defaults to (0,1)
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
#'    method = "ds", meta.data = list(width = 4))
#' plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1)
#' }
"plot_cond" <-function(obs,xmat,gxvalues,est,nc,breaks,showpoints,showlines,maintitle,ylim,angle=-45,density=20,col="black",jitter=NULL,xlab="Distance",ylab="Detection probability",...)
{
   selection <-xmat$detected[xmat$observer!=obs]==1
   selmat <- (xmat[xmat$observer==obs,])[selection,]
   shist<-hist(xmat$distance[xmat$observer!= obs & xmat$detected==1], breaks=breaks, plot = FALSE)
   mhist<-hist(xmat$distance[xmat$timesdetected== 2 & xmat$observer==obs], breaks=breaks, plot = FALSE)
   if(length(mhist$counts) < length(shist$counts)) 
     prop <- c(mhist$counts/shist$counts[1:length(mhist$counts)],rep(0, (length(shist$counts) - length(mhist$counts))))
   else 
     prop <- mhist$counts/shist$counts
   mhist$density<-prop
   mhist$equidist<-FALSE
   mhist$intensities<-mhist$density
   histline(mhist$density,breaks=breaks,lineonly=FALSE,xlab=xlab,ylab=ylab,ylim=ylim,
		     fill=TRUE,angle=angle,density=density,col=col,det.plot=TRUE,...) 
   if(showlines)
      lines(est$x,est$p,...)
   if(showpoints)
   {
	   ifelse(is.null(jitter),jitter.p<-1,jitter.p<-rnorm(length(gxvalues),1,jitter))
	   points(selmat$distance,gxvalues*jitter.p,...)
   }
   if(maintitle!="")
      title(paste(maintitle, "\nConditional detection probability - Observer=",obs ," | Observer = ",3-obs),...)
}
