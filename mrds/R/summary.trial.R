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
#' @S3method summary trial
#' @method summary trial
#' @aliases summary.trial
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects 
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.trial <-
function(object,se=TRUE,...)
{
#
# summmary.trial
#
# Provides a summary of parameters and estimates from the output of trial object
#
# Arguments:
#
# object     - object from ddf.trial
#
# Value: list of class summary.trial
#
# Uses: predict.trial, summary.ds, summary.trial.fi
#
model=object
avgp=function(model,pdot,...)
{return(pdot)}

  n=nrow(model$ds$ds$aux$ddfobj$xmat)
  ans=list(mr.summary=summary(model$mr,se=se,N=FALSE,fittedmodel=model),
           ds.summary=summary(model$ds,se=se, N=FALSE),
           Nhat=model$Nhat,AIC=model$criterion,
           average.p=n/model$Nhat)
#
# 26 Jan 06 jll; added code for se of average p
#
  if(se)
  {
     vcov=solvecov(model$hessian)$inv
     Nhatvar.list=DeltaMethod(model$par,NCovered,vcov,.001,model=model,group=TRUE)
     Nhatvar=Nhatvar.list$variance + sum((1-model$fitted)/model$fitted^2)
     cvN=sqrt(Nhatvar)/model$Nhat
     var.pbar.list=prob.se(model,avgp,vcov)
     covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
     var.pbar=ans$average.p^2*(cvN^2 + var.pbar.list$var/n^2
                                        -2*covar/(n*model$Nhat))
     ans$average.p.se=sqrt(var.pbar)
     ans$Nhat.se=sqrt(Nhatvar)
  }

  ans$mono<-model$ds$aux$mono
  ans$mono.strict<-model$ds$aux$mono.strict

  class(ans)="summary.trial"
  return(ans)
}
