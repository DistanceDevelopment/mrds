plot.det.tables=function(x,which=1:2,angle=-45,density=20,col1="black",col2="blue", new=TRUE,...)
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
	if(is.element(1,which)) 
	{
		if(new)dev.new()
		plot_seen(x$Observer1,"black","blue",...)
		legend("topright",legend=c("Detected by either observer","Detected by observer 1"),lty=1,lwd=3,col=c(col1,col2))
	}
	if(is.element(2,which)) 
	{
		if(new)dev.new()
		plot_seen(x$Observer2,"black","blue",...)
		legend("topright",legend=c("Detected by either observer","Detected by observer 2"),lty=1,lwd=3,col=c(col1,col2))
	}
	if(is.element(3,which)) 
	{
		if(new)dev.new()
		histline(x$Duplicates,breaks=breaks,lineonly=FALSE,xlab="Distance",ylab="Frequency",fill=TRUE,angle=angle,density=density,col=col1,...)
		legend("topright",legend=c("Seen by both observers"),lty=1,lwd=3,col=c(col1))
	}
	if(is.element(4,which)) 
	{
		if(new)dev.new()
		histline(x$Pooled,breaks=breaks,lineonly=FALSE,xlab="Distance",ylab="Frequency",fill=TRUE,angle=angle,density=density,col=col1,...)
		legend("topright",legend=c("Seen by either observer"),lty=1,lwd=3,col=c(col1))
	}
	if(is.element(5,which)) 
	{
		if(new)dev.new()
		plot_seen(x$Obs1_2,"black","blue",...)
		legend("topright",legend=c("Detected by observer 2","Detected by observer 1 | 2"),lty=1,lwd=3,col=c(col1,col2))
	}
	if(is.element(6,which)) 
	{
		if(new)dev.new()
		plot_seen(x$Obs2_1,"black","blue",...)
		legend("topright",legend=c("Detected by observer 1","Detected by observer 2 | 1"),lty=1,lwd=3,col=c(col1,col2))
	}
	return(NULL)
}