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
#' @S3method summary io
#' @method summary io
#' @aliases summary.io
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
summary.io <-
function(object,se=TRUE,...)
{
#
# summmary.io
#
# Provides a summary of parameters and estimates from the output of io object
#
# Arguments:
#
# object      - object from ddf.io
#
# Value: list of class summary.io.fi
#
# Uses: predict.io, summary.ds, summary.io.fi
#
model=object
avgp=function(model,pdot,...)
{return(pdot)}
  ddfobj=model$ds$ds$aux$ddfobj
  n=nrow(ddfobj$xmat)
  ans=list(mr.summary=summary(model$mr,se=se,N=FALSE,model,ddfobj),
           ds.summary=summary(model$ds,se=se, N=FALSE),
           Nhat=model$Nhat,AIC=model$criterion,
           average.p=n/model$Nhat,
           mono=model$ds$aux$mono,mono.strict=model$ds$aux$mono.strict)
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
     ans$cv=cvN
  }
  class(ans)="summary.io"
  return(ans)
}
