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
#' @S3method summary rem.fi
#' @method summary rem.fi
#' @aliases summary.rem.fi
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param N if TRUE, computes abundance in covered (sampled) region
#' @param fittedmodel full fitted model when called from \code{trial} or
#'   \code{io}
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects 
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.rem.fi <-
function(object,se=TRUE,N=TRUE,fittedmodel=NULL,...)
{
#
# summmary.trial.fi
#
# Provides a summary of parameters and estimates from the output of trial.fi object
#
# Arguments:
#
# object     - object from ddf.trial.fi
#
# Value: list of class summary.trial.fi
#
# Uses: predict.trial.fi
#
model=object
avgp=function(model,pdot,...)
{return(pdot)}

  newdat=model$data
  newdat <- newdat[newdat$distance <= model$meta.data$width &
        newdat$distance >= model$meta.data$left, ]
  n = length(newdat$distance)/2
  timesdetected <- newdat$detected[newdat$observer == 1] +
        newdat$detected[newdat$observer == 2]
  n1 = length(newdat$distance[newdat$observer == 1 & newdat$detected == 1])
  n2 = length(newdat$distance[newdat$observer == 2 & newdat$detected == 1])
  n3 = sum(timesdetected == 2)
  newdat$distance<-rep(0,length(newdat$distance))
  pred.at0=predict(model,newdat)$fitted
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
#
# Output conditional detection function values
#
  average.p0.1=sum(avgp(model,pred.at0)/pdot)

  ans=list(n=n,n1=n1,n2=n2,n3=n3, average.p0.1=average.p0.1/Nhat,cond.det.coef=coef(model),aic=model$criterion)
#
#   Output estimates relevant to abundance estimation.  Note: that these estimates are 
#   for the covered region only and are not expanded to some larger region that was sampled.
#
  if(N)
  {
     ans$average.p=n1/Nhat
     ans$Nhat=Nhat
  }
#
# Compute se for N and p's;  The logical N is set when summary.trial.fi is being
# called for method="trial.fi".  When method="trial" the call to summary.trial.fi from
# summary.trial sets N=FALSE.  So this is used to determine whether the pdots are to
# be computed from the trial.fi model or the trial model results.  Likewise when
# method="trial", fittedmodel contains the trial model object and is not NULL.
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
#    Compute se of p(0)
     if(is.null(fittedmodel))
        var.pbar.list=prob.se(model,avgp,vcov,fittedmodel=NULL)
     else
        var.pbar.list=prob.se(model,avgp,vcov,fittedmodel=fittedmodel)
     covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
     var.pbar=ans$average.p0.1^2*(cvN^2 + var.pbar.list$var/average.p0.1^2
                                     -2*covar/(average.p0.1*Nhat))
     ans$average.p0.1.se=sqrt(var.pbar)
  }

  ans$mono<-model$ds$aux$mono
  ans$mono.strict<-model$ds$aux$mono.strict

  class(ans)="summary.rem.fi"
  return(ans)
}
