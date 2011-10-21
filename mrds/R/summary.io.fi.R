#' Summary of distance detection function model object
#' 
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error.
#' 
#' The argument \code{N} is used to suppress computation of
#' abundance and average detection probability in calls to summarize the
#' \code{ds} and either \code{io.fi} or \code{trial.fi} for summaries of
#' \code{io} and \code{trial} objects respectively which are composed of a
#' \code{ds} model object and a mark-recapture model object. The corresponding
#' print function is called to print the summary results.
#' 
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param N if TRUE, computes abundance in covered (sampled) region
#' @param fittedmodel full fitted model when called from \code{trial} or
#'   \code{io}
#' @param ddfobj distance sampling object description
#' @param \dots unspecified and unused arguments for S3 consistency
#' @S3method summary io.fi
#' @method summary io.fi
#' @aliases summary.io.fi
#' @return list of extracted and summarized objects 
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.io.fi <- function(object,se=TRUE,N=TRUE,fittedmodel=NULL,ddfobj=NULL,...)
{
#
# summmary.io.fi
#
# Provides a summary of parameters and estimates from the output of io.fi object
#
# Arguments:
#
# object      - object from ddf.io.fi
#
# Value: list of class summary.io.fi
#
# Uses: predict.io.fi
#
model=object
avgp=function(model,pdot,...)
{return(pdot)}
avgp0=function(model,pdot,observer,...)
{
if(observer==1)
  return(pdot$p1)
else
  if(observer==2)
    return(pdot$p2)
  else
    return(pdot$fitted)
}

  newdat<-model$mr$data
  newdat$distance<-rep(0,length(newdat$distance))
  if(!is.null(ddfobj) && ddfobj$type=="gamma")
  {
	  key.scale <- scalevalue(ddfobj$scale$parameters,ddfobj$scale$dm)
	  key.shape <- scalevalue(ddfobj$shape$parameters,ddfobj$shape$dm)
	  newdat$distance=rep(apex.gamma(key.scale,key.shape),2)
  }  
  newdat$offsetvalue=0
  pred.at0=predict(model,newdat)
  if(is.null(fittedmodel))
  {
     Nhat=model$Nhat
     pdot=model$fitted
  }
  else
  {
    pdot=fittedmodel$fitted
    Nhat=fittedmodel$Nhat
  }
  n=length(newdat$distance)/2
  timesdetected <- newdat$detected[newdat$observer==1] + newdat$detected[newdat$observer==2]
  n1=length(newdat$distance[newdat$observer==1&newdat$detected==1])
  n2=length(newdat$distance[newdat$observer==2&newdat$detected==1])
  n3=sum(timesdetected==2)
  average.p0.1=sum(avgp0(model,pred.at0,observer=1)/pdot)
  average.p0.2=sum(avgp0(model,pred.at0,observer=2)/pdot)
  average.p0=sum(avgp0(model,pred.at0,observer=3)/pdot)

  ans=list(n=n,n1=n1,n2=n2,n3=n3,average.p0.1=average.p0.1/Nhat,average.p0.2=average.p0.2/Nhat,
           average.p0=average.p0/Nhat,cond.det.coef=coef(model),aic=model$criterion)
#
#   Output estimates relevant to abundance estimation.  Note: that these estimates are
#   for the surveyed strip only and are not expanded to some larger region that was sampled.
#
  if(N)
  {
     ans$average.p=n/Nhat
     ans$Nhat=Nhat
  }
#
# Compute se for N and p's;  The logical N is set when summary.io.fi is being
# called for method="io.fi".  When method="io" the call to summary.io.fi from
# summary.io sets N=FALSE.  So this is used to determine whether the pdots are to
# be computed from the io.fi model or the io model results.  Likewise when
# method="io", fittedmodel contains the io model object and is not NULL.
#
  if(se)
  {
     if(N | is.null(fittedmodel))
     {
        vcov=solvecov(model$hessian)$inv
        Nhatvar.list=DeltaMethod(model$par,NCovered,vcov,.001,model=model,group=TRUE)
        Nhatvar = Nhatvar.list$variance + sum((1-model$fitted)/model$fitted^2)
        cvN=sqrt(Nhatvar)/Nhat
        var.pbar.list=prob.se(model,avgp,vcov)
        covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
        var.pbar=ans$average.p^2*(cvN^2 + var.pbar.list$var/n^2
                                     -2*covar/(n*Nhat))
        ans$average.p.se=sqrt(var.pbar)
        ans$Nhat.se=sqrt(Nhatvar)
        ans$cv=cvN
     }
#    If there is a nested model, compute the Nhat from the fitted model and use
#    its variance and vector of partials.
     if(!is.null(fittedmodel))
     {
        Nhat=fittedmodel$Nhat
        vcov=solvecov(fittedmodel$hessian)$inv
        Nhatvar.list=DeltaMethod(fittedmodel$par,NCovered,vcov,.001,model=fittedmodel,group=TRUE)
        Nhatvar = Nhatvar.list$variance + sum((1-fittedmodel$fitted)/fittedmodel$fitted^2)
        cvN=sqrt(Nhatvar)/Nhat
     }
#    Compute se of p_1(0), p_2(0) and p_1(0)+p_2(0)-p_2(0)*p_1(0)
     var.pbar.list=prob.se(model,avgp0,vcov,observer=1,fittedmodel)
     covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
     var.pbar=ans$average.p0.1^2*(cvN^2 + var.pbar.list$var/average.p0.1^2
                                     -2*covar/(average.p0.1*Nhat))
     ans$average.p0.1.se=sqrt(var.pbar)
     var.pbar.list=prob.se(model,avgp0,vcov,observer=2,fittedmodel)
     covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
     var.pbar=ans$average.p0.2^2*(cvN^2 + var.pbar.list$var/average.p0.2^2
                                     -2*covar/(average.p0.2*Nhat))
     ans$average.p0.2.se=sqrt(var.pbar)
     var.pbar.list=prob.se(model,avgp0,vcov,observer=3,fittedmodel)
     covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
     var.pbar=ans$average.p0^2*(cvN^2 + var.pbar.list$var/average.p0^2
                                     -2*covar/(average.p0*Nhat))
     ans$average.p0.se=sqrt(var.pbar)

  }

  ans$mono<-model$ds$aux$mono
  ans$mono.strict<-model$ds$aux$mono.strict

  class(ans)="summary.io.fi"
  return(ans)
}
