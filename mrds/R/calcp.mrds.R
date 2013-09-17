#' Average detection probability calculation- DEPRECATED FUNCTION NO LONGER USED
#' 
#' Detection probability at a distance can vary depending on the other
#' ancillary covariates collected for an observation.  To plot a single
#' detection line, it is necessary to compute the expected probability value at
#' distance across the estimated proportions of the covariates in the
#' population.  This function does that calculation by observer, for duplicates
#' and cummulative for both observers.
#' 
#' @param dpformula valid formula for binary regression model detection
#'   function
#' @param dplink link function
#' @param dppars vector of detection function parameters, named in agreement
#'   with dpformula, and in same order (not sure if order is crucial)
#' @param dpdata data frame in format for ddf() mrds analysis
#' @param vname name of distance variable for x-axis (would be "distance" in
#'   usual ddf() data frame)
#' @param lower lower bound for plotting on x-axis
#' @param upper upper bound for plotting on x-axis
#' @param divisions number of divisions to use for numerical integration
#' @param type "line" for line transects, "point" for point transects (affects
#'   how integration is done)
#' @param objname name of variable denoting detection/sighting/object number
#'   (would be "object" in usual ddf() data frame)
#' @param obsname name of variable denoting observer/platform number (1 or 2)
#'   (would be "observer" in usual ddf() data frame)
#' @return \item{x}{distances at which probabilities are calculated}
#'   \item{p1}{probability for observer 1} \item{p2}{probability for observer
#'   2} \item{p3}{probability seen by both observers} \item{p.}{probability
#'   seen by at least one observer}
#' @author Who wrote this?
#' @references See Laake and Borchers (2004) chapter 6 in Advanced Distance
#'   Sampling by Buckland et al.
calcp.mrds<-function(dpformula,dplink,dppars,dpdata,vname,lower=0,upper,divisions=30,type=TRUE,objname="sighting",obsname="plat")
{
#  dpformula	valid formula for binary regression model detection function
#  dplink	link function
#  dppars	vector of detection function parameters, named in agreement with dpformula, and in same order (not sure if order is crucial)
#  dpdata	data frame in format for ddf() mrds analysis
#  vname	name of distance variable for x-axis (would be "distance" in usual ddf() data frame)
#  lower	lower bound for plotting on x-axis
#  upper	upper bound for plotting on x-axis
#  divisions	number of divisions to use for numerical integration
#  type		"line" for line transects, "point" for point transects (affects how integration is done)
#  objname	name of variable denoting detection/sighting/object number (would be "object" in usual ddf() data frame)
#  obsname	name of variable denoting observer/platform number (1 or 2) (would be "observer" in usual ddf() data frame)
# Define functions replicate.dat and closest (not used at present) 
	replicate.dat<-function (dat,vname,lower=0,upper,divisions=20)
#---------------------------------------------------------------------
# Function: replicate.dat
#
# Description: This function replicates rows of data in dat
#              divisions times, setting the vector of dat[,vname]s
#              for each replicate of a single sighting equal to a
#              sequence from lower to upper, at
#              equally spaced intervals. Also creates column of weigths for
#              summing (i.e. numerical integration).
#
#              The idea is that replicates and weights are used in numerical
#              integration over dat[,vname] for each sighting.
#
# Arguments:
#
# vname     = name of variable over which to integrate
# dat       = data frame
# divisions = number of intervals over which numerical integration
#             should be carried out
# lower     = lower limit for numerical integration
# upper     = maximum limit for numerical integration
#
#---------------------------------------------------------------------
	{
		nrep=divisions+1
		sumwt=rep(1,nrep)
		sumwt[1] = sumwt[nrep] = 0.5
		nsit2 = length(dat[[1]])
		d = list()
		for(i in names(dat)) {
			if(is.factor(dat[[i]])) {
				d[[i]]=factor(rep(dat[[i]],rep(nrep,nsit2)),levels=levels(dat[[i]]))
			} else {
				d[[i]]=rep(dat[[i]],rep(nrep,nsit2))
			}
		}
		d=as.data.frame(d)
		lower=lower
		upper=upper
		xvals=seq(lower,upper,length=nrep)
		d[,vname]=rep(xvals,nsit2)
		d$integration.weight=rep(sumwt,nsit2)
		d$index=rep((1:nrep),nsit2)
		return(list(dat=d,x=xvals))
	}
	
	closest=function(of,to,select="first")
#-------------------------------------------------------------------------------
# Returns index of closest element(s) in "of" to "to".
# "of" must be scalar or vector numeric
# "to" must be scalar
#-------------------------------------------------------------------------------
	{
		dif=abs(of-to)
		ind=which(dif==min(dif))
		if(select=="first") ind=min(which(dif==min(dif)))
		if(select=="last") ind=max(which(dif==min(dif)))
		return(ind)
	}
	
 type=ifelse(type,"line","point")
 if(type!="line" & type!="point") stop("type must be line or point")

 w=upper-lower
 binwidth=w/divisions
# dat<-dpdata[dpdata$plat!=0,all.vars(dpformula)]
 dat<-dpdata[!is.na(dpdata[,obsname]) & dpdata[,obsname]!=0,]
 replist=replicate.dat(dat,vname=vname,lower=lower,upper=upper,divisions=divisions)
 repdat=replist$dat
 x=replist$x
 repdat$p=p.det(dpformula,dplink,dppars,repdat)

 rn<-eval(parse(text=paste("repdat$",obsname,sep="")))
 rd1=repdat[rn==1,]
 rd2=repdat[rn==2,]
 rd3=rd.=rd1
 rd3$p=rd1$p*rd2$p
if(length(rd1$p)>0 & length(rd2$p)>0){
 rd.$p=rd1$p+rd2$p-rd3$p
}
else{
	if(length(rd1$p)>0){
		rd.$p<-rd1$p
	}
	else{
		if(length(rd2$p)>0){
			rd.$p<-rd2$p
		}
	}
}

 if(type=="line") {
   dp1=rd1$integration.weight*binwidth*rd1$p/w
   dp2=rd2$integration.weight*binwidth*rd2$p/w
   dp3=rd3$integration.weight*binwidth*rd3$p/w
   dp.=rd.$integration.weight*binwidth*rd.$p/w
 }
 if(type=="point") {
   dp1=rd1$integration.weight*binwidth*rd1$p*2*rd1[,vname]/w^2
   dp2=rd2$integration.weight*binwidth*rd2$p*2*rd2[,vname]/w^2
   dp3=rd3$integration.weight*binwidth*rd3$p*2*rd3[,vname]/w^2
   dp.=rd.$integration.weight*binwidth*rd.$p*2*rd.[,vname]/w^2
 }
# avg1=aggregate(dp1, by=list(sighting=rd1$sighting), sum)
# avg2=aggregate(dp2, by=list(sighting=rd2$sighting), sum)
# avg3=aggregate(dp3, by=list(sighting=rd3$sighting), sum)
 avg.=aggregate(dp., by=list(object=rd.[,objname]), sum)
# names(avg1)[2]=names(avg2)[2]=names(avg3)[2]=names(avg.)[2]="p"
 names(avg.)[2]="p"
 Pwt=(1/avg.$p)/(sum(1/avg.$p)) # only use best estimate of abundance of each kine (i.e. 1/P.) for weighting
 # put Pwt's back in replicate data:
 pwt=rep(Pwt,times=rep((divisions+1),length(avg.$p)))

 # (there is an issue above of whether to filter out those not seen by plat 1, plat 2 before weighting and averaging)
if(length(rd1$p)>0){
	p1<-aggregate(rd1$p*pwt,by=list(index=rd1$index),sum)
	names(p1)[2]<-"p"
} else{p1<-NULL}

if(length(rd2$p)>0){p2<-aggregate(rd2$p*pwt,by=list(index=rd2$index),sum)
	names(p2)[2]<-"p"
} else{p2<-NULL}

if(length(rd3$p)>0){ p3=aggregate(rd3$p*pwt,by=list(index=rd3$index),sum)
	names(p3)[2]<-"p"
} else{p3<-NULL}

if(length(rd.$p)>0){ p.=aggregate(rd.$p*pwt,by=list(index=rd.$index),sum)
	names(p.)[2]<-"p"
} else{p.<-NULL}
 return(list(x=x,p1=p1$p,p2=p2$p,p3=p3$p,p.=p.$p))
}
