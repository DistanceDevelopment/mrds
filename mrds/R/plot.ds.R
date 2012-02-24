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
#' @aliases plot.ds
#' @S3method plot ds
#' @method plot ds
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced. See
#'   \code{details.}
#' @param byvar name of variable to be used to color points - not currently
#'   implemented.
#' @param breaks user define breakpoints
#' @param nc number of equal-width bins for histogram
#' @param winht plot window height (not currently implemented)
#' @param winwd plot window width (not currently implemented)
#' @param jitter.v scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.v[j].  Where j=1,2 corresponds to
#'   observer j and j=3 corresponds to pooled/duplicate detections.
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param subset subset of data to plot
#' @param pl.col colours plotting colours for obs 1, obs 2 detections
#' @param bw.col grayscale plotting colours for obs 1, obs 2 detections
#' @param black.white logical variable; if TRUE plots are grayscale
#' @param pl.den shading density for plots of obs 1, obs 2 detections
#' @param pl.ang shading angle for plots of obs 1, obs 2 detections
#' @param main user-specfied plot title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (plot, hist, lines, points, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot
#' @examples
#' 
#' data(book.tee.data)
#' region<<-book.tee.data$book.tee.region
#' egdata<<-book.tee.data$book.tee.dataframe
#' samples<<-book.tee.data$book.tee.samples
#' obs<<-book.tee.data$book.tee.obs
#' xx=ddf(dsmodel = ~mcds(key = "hn", formula = ~sex), data = egdata[egdata$observer==1, ], method = "ds", meta.data = list(width = 4))
#' plot(xx,breaks=c(0,.5,1,2,3,4),showpoints=FALSE)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==0)
#' plot(xx,breaks=c(0,.5,1,2,3,4),subset=sex==1)
#' 
plot.ds <- function(x,which=c(2),byvar="",breaks=NULL,nc=NULL,winht=4,winwd=6,jitter.v=rep(0,3),showpoints=TRUE,subset=NULL,pl.col=c('black'),bw.col=c(grey(0)),black.white=FALSE,pl.den=rep(20,1),pl.ang=rep(-45,1),main=NULL,...)
{
# @usage \method{plot}{ds}(x,which=c(2),byvar="",breaks=NULL,nc=NULL,winht=4,winwd=6,jitter.v=rep(0,3),showpoints=TRUE,subset=NULL,pl.col=c('black'),bw.col=c(grey(0)),black.white=FALSE,pl.den=rep(20,1),pl.ang=rep(-45,1),...)
# plot.ds - detection function plot showing a scaled histogram of detections and then a line giving being the detection function
#           averaged over the estimated population levels of the covariate values, and one dot for each observation at its
#           estimated detection probability.
#
# Arguments:
# x                     - fitted model from ddf
# which                 - index to specify which plots should be produced:
#                             1 - data summary plot
#                             2 - data summary plot with overlaid average detection function and estimated detection probabilities
# byvar                 - name of variable to be used to color points - not currently implemented.
# breaks                - user define breakpoints
# nc                    - number of equal-width bins for histogram
# winht                 -	plot window height (not currently implemented)
# winwd	                - plot window width  (not currently implemented)
# jitter.v              - scaling option for plotting points.  Jitter is applied to points by multiplying the fitted value by a random
#                          draw from a normal distribution with mean 1 and sd jitter.v[j].  Where j=1,2 corresponds to observer j
#                          and j=3 corresponds to pooled/duplicate detections.
# pl.col                - colours plotting colours for obs 1, obs 2 detections
# bw.col                - grayscale plotting colours for obs 1, obs 2 detections
# black.white           - logical variable; if TRUE plots are grayscale
# pl.den                - shading density for plots of obs 1, obs 2 detections
# pl.ang                - shading angle for plots of obs 1, obs 2 detections
# showpoints            - logical variable; if TRUE plots predicted value for each observation
# subset                - subset of data to plot
# ...                   - other graphical parameters, passed to the plotting functions
#                          (plot, hist, lines, points, etc)
#
#  Uses: setcov, detfct, histline, test.breaks.
  model<-x
  show <- rep(FALSE, 8)
  show[which] <- TRUE
  lower<-0
  divisions<-25
  vname<-"distance"
  dpdata<-model$data
  dspars<-model$ds$par
#  dpdata$offsetvalue<-0

#est<-calcp.mrds(dpformula=dpformula,dplink=dplink,dppars=dppars,dpdata=dpdata,vname=vname,lower=lower,upper=upper,divisions=divisions,type="line",objname="object",obsname="observer")

  dspars<-model$ds$par
  dsmodel<-model$call$dsmodel[[2]][[2]]
  
  xlab<-"Distance"
  p1.name<-"Observer 1"
  p2.name<-"Observer 2"
  objname<-"object"
  obsname<-"observer"
  detname<-"detected"
  #nclass<-8

  dat<-dpdata

 # code from dpexplot:
  ltmodel<-model$ds
  width<-model$meta.data$width
  left<-model$meta.data$left
  ddfobj <- ltmodel$aux$ddfobj
  point=ltmodel$aux$point
  if(is.null(ltmodel$aux$int.range)){
    int.range=c(0,width)
  }else{
    int.range=ltmodel$aux$int.range
  }

  if(is.matrix(int.range)){
    max.range=as.vector(int.range[1,])
    int.range=int.range[2:dim(int.range)[1],]
    range.varies=TRUE
  }else{
    max.range=int.range
    normalize=FALSE
    range.varies=FALSE
  }

  if(range.varies&showpoints)
    errors("Point values can be misleading for g(x) when the range varies")

  if(!is.null(substitute(subset)))
    selected=eval(substitute(subset),ddfobj$xmat)
  else
    selected=rep(TRUE,nrow(ddfobj$xmat))

  if(all(!selected)){
    errors("Specified subset is empty.")
    return()
  }

  if(is.matrix(int.range)){
    int.range=int.range[selected,]
  }
  
  xmat<-ddfobj$xmat[selected,]
  if(!is.null(ddfobj$scale)){
    z<-ddfobj$scale$dm[selected,,drop=FALSE]
  }else{
    z<-matrix(1,nrow=1,ncol=1)
  }

  if(length(model$fitted)==1){
    pdot<-rep(model$fitted,sum(as.numeric(selected)))
  }else{
    pdot<-model$fitted[selected]
    Nhat <- sum(1/pdot)
  }

  zdim <- dim(z)[2]
  n <- length(xmat$distance)
  intercept.only <- ddfobj$intercept.only
  theta <- model$par

  if(!is.null(breaks)){
    nc<-length(breaks)-1
  }

  if(is.null(nc))
    nc<-round( sqrt(n),0)

#
#  Set logical hascov which if True means the detection function has covariates
#  other than distance and observer
#
  hascov <- FALSE
  if(!ddfobj$intercept.only){
    hascov <- TRUE
  }
#
#   Compute a grid for distance (xgrid), and covariates zgrid for
#   plotting of detection functions.  Also create intervals of distance (breaks)
#   for the chosen number of classes (nc).
#
  if(!hascov){
    xgrid <- (width/100)*(seq(0:100)-1)
    zgrid <- matrix(rep(z[1,],length(xgrid)), byrow=TRUE, ncol=sum(zdim))
  }

  if(is.null(breaks)){
    if(is.null(model$meta.data$binned)){
      binned<-FALSE
    }else{
      binned<-model$meta.data$binned
    }
    if(binned){
      breaks<-ltmodel$aux$breaks
      nc<-length(breaks)-1
    }else{
      breaks <-c(max(0,(max.range[1])), max.range[1] +((max.range[2]-max.range[1])/nc)*(1:nc))
      if(breaks[1]>left){
        breaks<-c(left,breaks)
        nc<-nc+1
      }
    }
  }
#
#   test breaks for validity and reset as needed
#
  breaks<-test.breaks(breaks,model$meta.data$left,width)
  nc<-length(breaks)-1

  #breaks<-seq(lower,upper,length=(nclass+1))
  lower<-min(breaks)
  upper<-max(breaks)
  dat<-dpdata[selected,]
  keep<-dat[,vname]>=lower & dat[,vname]<=upper
  h1<-hist(dat[,vname][keep],breaks=breaks,plot=FALSE)
  ymax<-max(h1$counts)

 #Set printing options for plots:
# By default  pl.col=c('black'),
#             bw.col=c(grey(0))

  # If greyscale plots are required use the following colours:
  if(black.white){
    byval1<-bw.col[1]
  }else{
  # If colour plots are required use the following:
    byval1=pl.col[1]
  }
 
  # Density of shaded lines - default is set all to 20
  denval1<-pl.den[1]
  # Angle of shaded lines - default is set all to -45
  angval1<-pl.ang[1]

  #Scaling for grouped data:
  if(normalize&!point){
    bindata<-function(x,r,breaks){
      return(hist(r[r>=x[1]&r<=x[2]],breaks=breaks,plot=FALSE)$counts)
    }
    sumit<-function(x,n,wt){
      return(sum(x/(wt*n)))
    }
    expected.counts<-apply(int.range,1,bindata,r=(0:1000)*width/1001,breaks=breaks)
    expected.counts<-apply(expected.counts,1,sumit,n=1001,wt=pdot)
  }else{
    if(!point){
      expected.counts<-((breaks[2:(nc+1)]-breaks[1:nc])*(Nhat/breaks[nc+1] ))
    }else{
	  expected.counts<--apply(matrix(c(breaks[2:(nc+1)]^2,breaks[1:nc]^2),ncol=2,nrow=nc),1,diff)*(Nhat/breaks[nc+1]^2 )
    }  
  }

  hist.obj<-hist(dat[,vname][keep], breaks=breaks, plot = FALSE)
  hist.obj$density<-hist.obj$counts/expected.counts
  hist.obj$density[expected.counts==0]=0
  freq<-hist.obj$density
  hist.obj$equidist<-FALSE
#------------------------------------------------------------------------


  # Data summary plot for observer 1
  if(show[1]){
    if(.Device=="windows"){
      win.graph()
    }
    histline(h1$counts,breaks=breaks,lineonly=FALSE,ylim=c(0,ymax),xlab=xlab,ylab="Frequency",fill=TRUE,angle=angval1,density=denval1,col=byval1,...)
    title(paste(p1.name,"detections"),cex.main=0.8)
  }

  # Detection function plot overlaid on historgram of observed 
  # distances - much of this code is taken from earlier plot.ds function.
  if(show[2]){
    if(.Device=="windows"){
      win.graph()
    }

    # Primary detection function
    gxvalues <- detfct(xmat$distance,ddfobj,select=selected,width=width)
    ymax <- max(hist.obj$density,max(gxvalues))
    histline(hist.obj$density,breaks=breaks,lineonly=FALSE,
             xlab=xlab,ylab="Detection probability",ylim=c(0,ymax),
             fill=TRUE,angle=angval1,density=denval1,col=byval1,
             det.plot=TRUE,...)

    if(hascov){
      finebr <- (width/100)*(0:100)
      xgrid <- NULL
      linevalues <- NULL
      newdat <- xmat
      for(i in 1:(length(finebr)-1)){
        x <- (finebr[i]+finebr[i+1])/2
        xgrid <- c(xgrid,x)
        newdat$distance <- rep(x,nrow(newdat))

        detfct.values <-detfct(newdat$distance,ddfobj,select=selected,
                               width=width)
        if(!normalize&range.varies){
          detfct.values[x<int.range[,1]|x>int.range[,2]]=0
        }
        linevalues<-c(linevalues,sum(detfct.values/pdot)/sum(1/pdot))
      }
    }else{	  
      if(!is.null(ddfobj$scale)){
        ddfobj$scale$dm <- ddfobj$scale$dm[rep(1,length(xgrid)),]
      }
	   if(!is.null(ddfobj$shape)){
        ddfobj$shape$dm <- ddfobj$shape$dm[rep(1,length(xgrid)),]
      }
	   linevalues <- detfct(xgrid,ddfobj,width=width)
    }
    lines(xgrid,linevalues,col=byval1,...)
    if(showpoints){
      jitter.p <- rnorm(length(gxvalues),1,jitter.v[1])
      points(xmat$distance,gxvalues*jitter.p,col=byval1,...)
    }

    # use the user-supplied main= ...
    if(is.null(main)){
       title("Detection function plot",cex.main=0.8)
    }else{
       title(main,cex.main=0.8)
    }
    invisible()
  }

}
