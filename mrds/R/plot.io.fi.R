

#' Plot fit of detection functions and histograms of data from distance
#' sampling model
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
#' @aliases plot.io.fi
#' @method plot io.fi
#' @S3method plot io.fi
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
#' @param showlines logical variable; if TRUE a line representing the average
#'   detection probability is plotted
#' @param pl.col colours plotting colours for obs 1, obs 2 detections
#' @param bw.col grayscale plotting colours for obs 1, obs 2 detections
#' @param black.white logical variable; if TRUE plots are grayscale
#' @param pl.den shading density for plots of obs 1, obs 2 detections
#' @param pl.ang shading angle for plots of obs 1, obs 2 detections
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
plot.io.fi <- function(x,which=c(1,2,7,8),byvar="",breaks=NULL,nc=NULL,winht=4,winwd=6,jitter.v=rep(0,3),showpoints=TRUE,showlines=TRUE,pl.col=c('black','blue','green','red'),bw.col=c(grey(0),grey(0.4),grey(0.8),grey(0.6)),black.white=FALSE,pl.den=rep(20,4),pl.ang=rep(-45,4),...){
# plot.io.fi - Provides plots of fitted functions for a io.fi object
#
# Arguments:
# x                     - fitted model from ddf
# which                 - index to specify which plots should be produced:
#                             1 - data summary plot for observer 1
#                             2 - data summary plot for observer 2
#                             3 - observer 1 detection function plot with MR model fit and data
#                             4 - observer 2 detection function plot with MR model fit and data
#                             5 - duplicate detection function plot with MR model fit and data
#                             6 - combined detection function plot with MR model fit and data
#                             7 - conditional MR detection function plot obs 1 given obs 2 giving proportion of duplicates with fitted MR model
#                                   averaged over population covariate values and dots for each estimated det prob.
#                             8 - conditional MR detection function plot obs 2 given obs 1 giving proportion of duplicates with fitted MR model
#                                   averaged over population covariate values and dots for each estimated det prob.
# byvar                 - name of variable to be used to color points - not currently implemented.
# breaks                - user define breakpoints
# nc                    - number of equal-width bins for histogram
# winht                 -	plot window height (not currently implemented)
# winwd	                - plot window width  (not currently implemented)
# jitter.v              - scaling option for plotting points.  Jitter is applied to points by multiplying the fitted value by a random
#                          draw from a normal distribution with mean 1 and sd jitter.v[j].  Where j=1,2 corresponds to observer j
#                          and j=3 corresponds to pooled/duplicate detections.
# showpoints            - logical variable; if TRUE plots predicted value for each observation
# showlines             - logical variable; if TRUE a line representing the average detection probability is plotted
# pl.col                - colours plotting colours for obs 1, obs 2, pooled and duplicate detections
# bw.col                - grayscale plotting colours for obs 1, obs 2, pooled and duplicate detections
# black.white           - logical variable; if TRUE plots are grayscale
# pl.den                - shading density for plots of obs 1, obs 2, pooled and duplicate detections
# pl.ang                - shading angle for plots of obs 1, obs 2, pooled and duplicate detections
# ...                   - other graphical parameters, passed to the plotting functions
#                          (plot, hist, lines, points, etc)
#
#  Uses: calcp.mrds, detfct, histline, average.line
 model<-x
 show <- rep(FALSE, 8)
 show[which] <- TRUE

  dpformula<-as.formula(model$mr$formula)
  dplink<-model$mr$family$link
#  upper<-model$call$meta.data$width
  upper<-model$meta.data$width
  dppars<-model$mr$coefficients
  regs<-all.vars(dpformula)
  lower<-0
  divisions<-25
  vname<-"distance"
#  dpdata<-model$mr$data - for trial.fi method this only contains observations from obs1
#  dpdata<-model$data
   dpdata<-model$mr$data
  dpdata$offsetvalue<-0


   xmat<- model$mr$data
   xmat$offsetvalue=0

   width <- model$meta.data$width
   left <- model$meta.data$left
   lower<-left
   
#dpdata$sizefactor=as.factor(ddf.all.dat$sizefactor) # make sure size is a factor
est<-calcp.mrds(dpformula=dpformula,dplink=dplink,dppars=dppars,dpdata=dpdata,vname=vname,lower=lower,upper=upper,divisions=divisions,type="line",objname="object",obsname="observer")
xlab<-"Distance"
p1.name<-"Observer 1"
p2.name<-"Observer 2"
objname<-"object"
obsname<-"observer"
detname<-"detected"
#nclass<-8


       if(is.null(nc))
    nc<-round(sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
                   length(xmat$distance[xmat$observer==2&xmat$detected==1]),length(xmat$distance[xmat$observer==1&xmat$timesdetected==2]) )),0)
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

# code from plot.mrds

 keep1=dpdata[,obsname]==1 & dpdata[,detname]==1
 p1=p.det(dpformula,dplink,dppars,dpdata[keep1,])
# r1=dpdata$r[keep1]
 r1=dpdata[keep1,vname]
 keep2=dpdata[,obsname]==2 & dpdata[,detname]==1
 p2=p.det(dpformula,dplink,dppars,dpdata[keep2,])
 r2=dpdata[keep2,vname]
 p1.=p.det(dpformula,dplink,dppars,dpdata[dpdata[,obsname]==1,])
 p2.=p.det(dpformula,dplink,dppars,dpdata[dpdata[,obsname]==2,])
 p3<-p1.*p2.

if(length(p1.)>0 & length(p2.)>0){
 p.=p1.+p2.-p1.*p2.
}
else{
	if(length(p1.)>0){
		p.=p1.
	}
	else{
		if(length(p2.)>0){
			p.=p2.
		}
	}
}
# r.=dpdata$r[dpdata[,obsname]==1]
 keep.=dpdata[,obsname]==1 # used in plot.mrds
 r.=dpdata[keep.,vname] #- used in plot.mrds

# get "byvar" variables
 if(byvar=="") {
   byval1=byval2=byval.=1
 } else {
   regs=all.vars(dpformula)
   if(!is.element(byvar,regs)) {
     byval=1
     warning(paste("Invalid byvar: byvar is",byvar," but variables are:",paste(regs,collapse=", ")))
   }else {
     byval1=as.numeric(dpdata[keep1,byvar])
     byval2=as.numeric(dpdata[keep2,byvar])
     byval.=as.numeric(dpdata[dpdata[,obsname]==1,byvar])
     if(min(byval1)<1) byval1=byval1+(1-min(byval1))
     if(min(byval2)<1) byval2=byval2+(1-min(byval2))
     if(min(byval.)<1) byval.=byval.+(1-min(byval.))
   }
 }
 dat=dpdata

 # code from dpexplot:
 s1=(dat[,obsname]==1 & dat[,detname]==1)
 s2=(dat[,obsname]==2 & dat[,detname]==1)
 sn1=dat[,objname][s1]
 sn2=dat[,objname][s2]
 snd=sn1[which(is.element(sn1,sn2))]

# breaks=seq(lower,upper,length=(nclass+1))

 keep=dat[,vname]>=lower & dat[,vname]<=upper
 h1=hist(dat[,vname][keep & is.element(dat[,objname],sn1) & dat[,obsname]==1],breaks=breaks,plot=FALSE)
 h2=hist(dat[,vname][keep & is.element(dat[,objname],sn2) & dat[,obsname]==2],breaks=breaks,plot=FALSE)
 hd1=hist(dat[,vname][keep & is.element(dat[,objname],snd) & dat[,obsname]==1],breaks=breaks,plot=FALSE)
 hd2=hist(dat[,vname][keep & is.element(dat[,objname],snd) & dat[,obsname]==2],breaks=breaks,plot=FALSE)

 ymax=max(h1$counts,h2$counts,hd1$counts,hd2$counts)
 
#Set printing options for plots:
# By default  pl.col=c('black','blue','green','red'),
#             bw.col=c(grey(0),grey(0.4),grey(0.8),grey(0.6))
#If greyscale plots are required use the following colours:

if(black.white){
 byval1=bw.col[1]
 byval2=bw.col[2]
 byval3=bw.col[3]
 byval4=bw.col[4]
}
#If colour plots are required use the following:
else{
 byval1=pl.col[1]
 byval2=pl.col[2]
 byval3=pl.col[3]
 byval4=pl.col[4]
}
#Density of shaded lines - default is set all to 20
denval1<-pl.den[1]
denval2<-pl.den[2]
denval3<-pl.den[3]
denval4<-pl.den[4]

#Angle of shaded lines - default is set all to -45
angval1<-pl.ang[1]
angval2<-pl.ang[2]
angval3<-pl.ang[3]
angval4<-pl.ang[4]

 #Data summary plot for observer 1
if(show[1]){
 if (.Device=="windows") win.graph()
 histline(hd1$counts,breaks=breaks,lineonly=FALSE,ylim=c(0,ymax),xlab=xlab,ylab="Frequency",fill=TRUE,angle=angval1,density=denval1,col=byval1,...)
 histline(h1$counts,breaks,lineonly=TRUE,col=byval2,...)
 title(paste(p1.name,"detections"),cex.main=0.8)
}
 #Data summary plot for observer 2
if(show[2]){
 if (.Device=="windows") win.graph()
 histline(hd2$counts,breaks=breaks,lineonly=FALSE,ylim=c(0,ymax),xlab=xlab,ylab="Frequency",fill=TRUE,angle=angval2,density=denval2,col=byval2,...)
 histline(h2$counts,breaks,lineonly=TRUE,col=byval1,...)
 title(paste(p2.name,"detections"),cex.main=0.8)
}
 dp1=hd1$count/h1$count
 dp2=hd2$count/h2$count

#Observer 1 detection function plot:
#observer 1 detections with MR model fit and data – can use to see if full indep model might be better
if(show[3]){
  if (.Device=="windows") win.graph()
  #Where is the code for this? use average line code - from plot.uncond?
#  plot.uncond - plots unconditional detection function for observer=obs observations
#              overlays histrogram, average detection function and values for individual observations

#  xmatmodel$data
  xmat<-model$mr$data
  xmat$offsetvalue<-0
  obs<-1
  gxvalues<-p1
#  nc<-nclass
  finebr<-seq(0,max(breaks),length=50)
  keep=dat[,vname]>=lower & dat[,vname]<=upper

  n <- length(xmat$distance[keep & xmat$observer == obs & xmat$detected==1])
  selmat <- xmat[keep & xmat$observer== obs,]
  det.detected <- xmat$detected[keep & xmat$observer== obs]==1
  
  hist.obj<-hist(selmat$distance[det.detected], breaks=breaks, plot = FALSE)
  expected.counts=((breaks[2:(nc+1)]-breaks[1:nc])*(model$Nhat/breaks[nc+1] ))
  hist.obj$density<-hist.obj$counts/(expected.counts)
  hist.obj$intensities<-hist.obj$density
  freq<-hist.obj$density
  hist.obj$equidist<-FALSE
  line=average.line(finebr,obs,model)
  linevalues <- line$values
  xgrid <- line$xgrid
  maxl <- max(freq, gxvalues)

#  ylim<-c(0,max(ylim,hist.obj$density)) # *** DLB change ***
  ylim<-c(0,max(linevalues,maxl)) # *** JRBB change ***
#  plot(hist.obj,freq=FALSE,main="",ylab="Detection probability",xlab="Distance",xlim = c(0, model$meta.data$width),ylim=ylim,fill=TRUE,angle=-45,density=20,col=byval1)
  histline(hist.obj$density,breaks=breaks,lineonly=FALSE,ylab="Detection probability",xlab="Distance",ylim=ylim,fill=TRUE,angle=angval1,density=denval1,col=byval1,det.plot=TRUE,...)
  if(showlines){lines(xgrid,linevalues,col=byval1,...)}
  #  points(selmat$distance[det.detected],gxvalues)
  if(showpoints){
  jitter.p<-rnorm(length(gxvalues),1,jitter.v[1])
  points(selmat$distance[det.detected],gxvalues*jitter.p,col=byval1,...)
  }
        title(paste("\nObserver",obs, " detections"),cex.main=0.8)
}

#Observer 2 detection function plot:
if(show[4]){
  if (.Device=="windows") win.graph()
  #Where is the code for this? use average line code - from plot.uncond?
#  plot.uncond - plots unconditional detection function for observer=obs observations
#              overlays histrogram, average detection function and values for individual observations

#  xmat<-model$data
  xmat<-model$mr$data
  xmat$offsetvalue<-0
  obs<-2
  gxvalues<-p2
 # nc<-nclass
  finebr<-seq(0,max(breaks),length=50)
  keep=dat[,vname]>=lower & dat[,vname]<=upper

  n <- length(xmat$distance[keep & xmat$observer == obs & xmat$detected==1])
  selmat <- xmat[keep & xmat$observer== obs,]
  det.detected <- xmat$detected[keep & xmat$observer== obs]==1

  hist.obj<-hist(selmat$distance[det.detected], breaks=breaks, plot = FALSE)
  expected.counts=((breaks[2:(nc+1)]-breaks[1:nc])*(model$Nhat/breaks[nc+1] ))
  hist.obj$density<-hist.obj$counts/(expected.counts)
  hist.obj$intensities<-hist.obj$density
  freq<-hist.obj$density
  hist.obj$equidist<-FALSE
  line=average.line(finebr,obs,model)
  linevalues <- line$values
  xgrid <- line$xgrid
  maxl <- max(freq, gxvalues)

#  ylim<-c(0,max(ylim,hist.obj$density)) # *** DLB change ***
  ylim<-c(0,max(linevalues,maxl)) # *** JRBB change ***
#  plot(hist.obj,freq=FALSE,main="",ylab="Detection probability",xlab="Distance",xlim = c(0, model$meta.data$width),ylim=ylim,fill=TRUE,angle=-45,density=20,col=byval1)
  histline(hist.obj$density,breaks=breaks,lineonly=FALSE,ylab="Detection probability",xlab="Distance",ylim=ylim,fill=TRUE,angle=angval2,density=denval2,col=byval2,det.plot=TRUE,...)
  if(showlines){lines(xgrid,linevalues,col=byval2,...)}
  #  points(selmat$distance[det.detected],gxvalues)
  if(showpoints){
  jitter.p<-rnorm(length(gxvalues),1,jitter.v[2])
  points(selmat$distance[det.detected],gxvalues*jitter.p,col=byval2,...)
  }
        title(paste("\nObserver",obs, " detections"),cex.main=0.8)
}

#5.duplicate detection function plot (obs=4)
if(show[5]){
  if (.Device=="windows") win.graph()
  #Where is the code for this? use average line code - from plot.uncond?
#  plot.uncond - plots unconditional detection function for observer=obs observations
#              overlays histrogram, average detection function and values for individual observations

  xmat<-model$data
  xmat.mr<-model$mr$data
  obs<-4
  gxvalues<-p3
 # nc<-nclass
  finebr<-seq(0,max(breaks),length=50)
  keep=dat[,vname]>=lower & dat[,vname]<=upper

  n <- length(xmat$distance[xmat$observer==1 & xmat.mr$timesdetected==2])
  selmat <- xmat.mr[xmat.mr$observer==1,]
  det.detected <- selmat$timesdetected==2

  hist.obj<-hist(selmat$distance[det.detected], breaks=breaks, plot = FALSE)
  expected.counts=((breaks[2:(nc+1)]-breaks[1:nc])*(model$Nhat/breaks[nc+1] ))
  hist.obj$density<-hist.obj$counts/(expected.counts)
  hist.obj$intensities<-hist.obj$density
  freq<-hist.obj$density
  hist.obj$equidist<-FALSE
  line=average.line(finebr,obs,model)
  linevalues <- line$values
  xgrid <- line$xgrid
  maxl <- max(freq, gxvalues)

#  ylim<-c(0,max(ylim,hist.obj$density)) # *** DLB change ***
  ylim<-c(0,max(linevalues,maxl)) # *** JRBB change ***
#  plot(hist.obj,freq=FALSE,main="",ylab="Detection probability",xlab="Distance",xlim = c(0, model$meta.data$width),ylim=ylim,fill=TRUE,angle=-45,density=20,col=byval1)
  histline(hist.obj$density,breaks=breaks,lineonly=FALSE,ylab="Detection probability",xlab="Distance",ylim=ylim,fill=TRUE,angle=angval4,density=denval4,col=byval4,det.plot=TRUE,...)
  if(showlines){lines(xgrid,linevalues,col=byval4,...)}
#   points(selmat$distance[det.detected],gxvalues)
  if(showpoints){
    jitter.p<-rnorm(length(gxvalues),1,jitter.v[3])
    points(selmat$distance,gxvalues*jitter.p,col=byval4,...)
  }
        title(paste("\n Duplicate detections"),cex.main=0.8)
}


#6.pooled detection function plot
if(show[6]){
  if (.Device=="windows") win.graph()
  #Where is the code for this? use average line code - from plot.uncond?
#  plot.uncond - plots unconditional detection function for observer=obs observations
#              overlays histrogram, average detection function and values for individual observations

#  xmat<-model$data
  xmat<-model$mr$data
  xmat$offsetvalue<-0
  obs<-3
  gxvalues<-p.
#  nc<-nclass
  finebr<-seq(0,max(breaks),length=50)
  keep=dat[,vname]>=lower & dat[,vname]<=upper

  n <- length(xmat$distance[xmat$observer==1])
  selmat <- xmat[xmat$observer==1,]
  det.detected<- selmat$observer==1

  hist.obj<-hist(selmat$distance[det.detected], breaks=breaks, plot = FALSE)
  expected.counts=((breaks[2:(nc+1)]-breaks[1:nc])*(model$Nhat/breaks[nc+1] ))
  hist.obj$density<-hist.obj$counts/(expected.counts)
  hist.obj$intensities<-hist.obj$density
  freq<-hist.obj$density
  hist.obj$equidist<-FALSE
  line=average.line(finebr,obs,model)
  linevalues <- line$values
  xgrid <- line$xgrid
  maxl <- max(freq, gxvalues)

#  ylim<-c(0,max(ylim,hist.obj$density)) # *** DLB change ***
  ylim<-c(0,max(linevalues,maxl)) # *** JRBB change ***
  histline(hist.obj$density,breaks=breaks,lineonly=FALSE,ylab="Detection probability",xlab="Distance",ylim=ylim,fill=TRUE,angle=angval3,density=denval3,col=byval3,det.plot=TRUE,...)
  if(showlines){lines(xgrid,linevalues,col=byval3,...)}
#    points(selmat$distance[det.detected],gxvalues)
  if(showpoints){
   xmat.mr<-model$mr$data
  #Seen by obs1 only:
  obs1.ind<-xmat$object[which(xmat$observer==1 & xmat$detected==1 & xmat.mr$timesdetected==1)]
  #Seen by obs2 only:
  obs2.ind<-xmat$object[which(xmat$observer==2 & xmat$detected==1 & xmat.mr$timesdetected==1)]
  #Seen by both obs:
  obsb.ind<-xmat$object[which(xmat$observer==1 & xmat$detected==1 & xmat.mr$timesdetected==2)]
    jitter.pb<-rnorm(length(gxvalues[obsb.ind]),1,jitter.v[3])
    jitter.p1<-rnorm(length(gxvalues[obs1.ind]),1,jitter.v[1])
    jitter.p2<-rnorm(length(gxvalues[obs2.ind]),1,jitter.v[2])
    points(selmat$distance[obsb.ind],gxvalues[obsb.ind]*jitter.pb,col=byval4,...)
    points(selmat$distance[obs1.ind],gxvalues[obs1.ind]*jitter.p1,col=byval1,...)
    points(selmat$distance[obs2.ind],gxvalues[obs2.ind]*jitter.p2,col=byval2,...)
}
#    points(selmat$distance,gxvalues,col=byval1,pch=19)
        title(paste("\nPooled detections"),cex.main=0.8)
}

#7.conditional detection function plot obs 1 | obs 2
if(show[7]){
  if (.Device=="windows") win.graph()
  if(length(p1.)>0 & length(p2.)>0){
    histline(dp2,breaks=breaks,lineonly=FALSE,ylim=c(0,1),xlab=xlab,ylab="Duplicate Proportion",fill=TRUE,angle=angval1,density=denval1,col=byval1,det.plot=TRUE,...)
#    title(paste("Proportion of ",p2.name,"detections seen by",p1.name),cex.main=0.75)
    title("Conditional MR detection plot Obs 1 | Obs 2",cex.main=0.8)
    if(showlines){lines(est$x,est$p1,type="l",ylim=c(0,1),col=byval1,...)}
    if(showpoints){
    jitter.p<-rnorm(length(p1),1,jitter.v[1])
    points(r1,p1*jitter.p,col=byval1,...)
    }
  }
}

#8.conditional detection function plot obs 2 | obs 1
if(show[8]){
  if (.Device=="windows") win.graph()
  if(length(p1.)>0 & length(p2.)>0){
    histline(dp1,breaks=breaks,lineonly=FALSE,ylim=c(0,1),xlab=xlab,ylab="Duplicate Proportion",fill=TRUE,angle=angval2,density=denval2,col=byval2,det.plot=TRUE,...)
#    title(paste("Proportion of ",p1.name,"detections seen by",p2.name),cex.main=0.75)
    title("Conditional MR detection plot Obs 2 | Obs 1",cex.main=0.8)
    if(showlines){ lines(est$x,est$p2,type="l",ylim=c(0,1),col=byval2,...)}
#    gxvalues <-p2[xmat$detected[xmat$observer==1]==1]
#    plot.cond(model,2,xmat,gxvalues,nc,finebr,breaks,showpoints,showlines,maintitle,ylim,...)
    if(showpoints){
    jitter.p<-rnorm(length(p2),1,jitter.v[2])
    points(r2,p2*jitter.p,col=byval2,...)
    }
  }
}

}