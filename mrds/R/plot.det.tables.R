#' Observation detection tables
#' 
#' Creates a series of tables for dual observer data that shows the number
#' missed and detected for each observer within defined distance classes.
#' 
#' @aliases plot.det.tables
#' @method plot det.tables
#' @S3method plot det.tables
#' @export
#' @param x object of class det.tables
#' @param which items in x to plot; vector with values in 1:6
#' #' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col1 plotting colour for specified universe of detections (col1= Observer1; col2= Observer 2 within Observer 1 subset) 
#' @param col2 plotting colour for those detected
#' @param new if TRUE new plotting window for each plot
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (plot, hist, lines, points, etc)
#' @return NULL
#' @author Jeff Laake
#' @examples
#' \donttest{
#' data(book.tee.data)
#' region<<-book.tee.data$book.tee.region
#' egdata<<-book.tee.data$book.tee.dataframe
#' samples<<-book.tee.data$book.tee.samples
#' obs<<-book.tee.data$book.tee.obs
#' xx=ddf(mrmodel=~glm(formula=~distance*observer),dsmodel = ~mcds(key = "hn", formula = ~sex), 
#'         data = egdata, method = "io", meta.data = list(width = 4))
#' tabs=det.tables(xx,breaks=c(0,.5,1,2,3,4))
#' par(mfrow=c(2,3))
#' plot(tabs,which=1:6,new=FALSE)
#' }
plot.det.tables=function(x,which=1:6,angle=-45,density=20,col1="black",col2="blue", new=TRUE,...)
{
	plot_seen=function(x,col1,col2,...)
	{
		missed=x[,"Missed"]
		detected=x[,"Detected"]
		ymax=max(missed+detected)
    	histline(detected,breaks=breaks,lineonly=FALSE,ylim=c(0,ymax),xlab="Distance",ylab="Frequency",fill=TRUE,angle=angle,density=density,col=col2,...)
	    histline(missed+detected,breaks,lineonly=TRUE,col=col1,...)
	}	
	breaks=x$breaks
	if(is.element(1,which)&!is.null(x$Observer1)) 
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_seen(x$Observer1,col1,col2,...)
		legend("topright",legend=c("Detected by either observer","Detected by observer 1"),lty=1,lwd=3,col=c(col1,col2))
	}
	if(is.element(2,which)&!is.null(x$Observer2)) 
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_seen(x$Observer2,col1,col2,...)
		legend("topright",legend=c("Detected by either observer","Detected by observer 2"),lty=1,lwd=3,col=c(col1,col2))
	}
	if(is.element(3,which)&!is.null(x$Duplicates)) 
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		histline(x$Duplicates,breaks=breaks,lineonly=FALSE,xlab="Distance",ylab="Frequency",fill=TRUE,angle=angle,density=density,col=col1,...)
		legend("topright",legend=c("Seen by both observers"),lty=1,lwd=3,col=c(col1))
	}
	if(is.element(4,which)&!is.null(x$Pooled))
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		histline(x$Pooled,breaks=breaks,lineonly=FALSE,xlab="Distance",ylab="Frequency",fill=TRUE,angle=angle,density=density,col=col1,...)
		legend("topright",legend=c("Seen by either observer"),lty=1,lwd=3,col=c(col1))
	}
	if(is.element(5,which)&!is.null(x$Obs1_2)) 
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_seen(x$Obs1_2,col1,col2,...)
		legend("topright",legend=c("Detected by observer 2","Detected by observer 1 | 2"),lty=1,lwd=3,col=c(col1,col2))
	}
	if(is.element(6,which)&!is.null(x$Obs2_1)) 
	{
		if(new& .Platform$GUI=="Rgui")dev.new()
		plot_seen(x$Obs2_1,col1,col2,...)
		legend("topright",legend=c("Detected by observer 1","Detected by observer 2 | 1"),lty=1,lwd=3,col=c(col1,col2))
	}
	return(NULL)
}
