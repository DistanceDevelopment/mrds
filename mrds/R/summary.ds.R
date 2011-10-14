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
#' @S3method summary ds
#' @method summary ds
#' @aliases summary.ds
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param N if TRUE, computes abundance in covered (sampled) region
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.ds <-
function(object,se=TRUE,N=TRUE,...)
{
#
# summmary.ds
#
# Provides a summary of parameters and estimates from a ds model object
#
# Arguments:
#
#  object  a ds model object
#  se      if TRUE, computes standard error of abundance in covered region
#  N       if TRUE, computes abundance in covered (sampled) region
#
# Value: Null
#
# Uses: predict.ds (via predict), DeltaMethod, coef.ds (via coef)
#
# dlm 31-Aug-05	Started re-working this to return an object.
# dlm 03-Sep-05	Will now return an object of type summary.ds, but will
# 		still behave as it used to. See print.summary.ds for details.
#
# at present f(0) code is commented out; need to understand this further
#F0=function(model,pdot,...)
#{return(1/(pdot*model$meta.data$width))}
model=object
avgp=function(model,pdot,...)
{return(pdot)}

  ans <- list()

# was monotonicity enforced?
  ans$mono<-model$ds$aux$mono
# strict monotonicity?
  ans$mono.strict<-model$ds$aux$mono.strict

# Number of observations
  ans$n <- length(model$ds$aux$ddfobj$xmat$distance)
# Average detection prob for mcds
  if(is.null(model$fitted))
    pdot=predict(model,esw=FALSE)$fitted
  else
    pdot=model$fitted

# Set the key function type
  ans$key <- model$ds$aux$ddfobj$type

# Parameter estimates and se for detection function
# se is included as part of the objects, see coef.ds
# for details
  coeff <- coef(model)

#  Scale Coefficient
  ans$coeff$key.scale <- coeff$scale

# Hazard shape parameter
  if(ans$key%in%c("gamma","hr")){
    ans$coeff$key.shape <- coeff$exponent
  }
  
# Adjustment term parameter(s)
# See coef.ds() on how this is returned
# This is a vector remember, so if you are using this 
# you need to take that into account.
  if(!is.null(coeff$adjustment)){
    ans$coeff$adj.order <- model$adj.order
    ans$coeff$adj.parm <- coeff$adjustment
  }
  

# AIC
  ans$aic <- model$criterion

# Truncation distances
# right width
  ans$width <- model$meta.data$width
# left width
  ans$left <- model$meta.data$left
  
  ans$average.p=ans$n/model$Nhat
#  fzero=sum(F0(model,pdot)/pdot)
#  ans$average.f0=fzero/model$Nhat
#
# 26 Jan 06 jll; added code for se of average p and f(0)
#
  if(se)
  {
     vcov=solvecov(model$hessian)$inv
     Nhatvar.list=DeltaMethod(model$par,NCovered,vcov,0.001,model=model,group=TRUE)
     Nhatvar = Nhatvar.list$variance + sum((1-model$fitted)/model$fitted^2)
     cvN=sqrt(Nhatvar)/model$Nhat
#     var.fzero.list=prob.se(model,F0,vcov)
#     covar=t(Nhatvar.list$partial)%*%vcov%*%var.fzero.list$partial+var.fzero.list$covar
#     var.fzero=ans$average.f0^2*(cvN^2 + var.fzero.list$var/fzero^2
#                                        -2*covar/(fzero*model$Nhat))
     var.pbar.list=prob.se(model,avgp,vcov)
     covar=t(Nhatvar.list$partial)%*%vcov%*%var.pbar.list$partial+var.pbar.list$covar
     var.pbar=ans$average.p^2*(cvN^2 + var.pbar.list$var/ans$n^2
                                        -2*covar/(ans$n*model$Nhat))
  }
  if(N)
  {
    ans$Nhat <- model$Nhat
    if(se)
       ans$Nhat.se <- sqrt(Nhatvar)
  }
  if(se)
  {
#    ans$average.f0.se=sqrt(var.fzero)
    ans$average.p.se=sqrt(var.pbar)
  }
  class(ans) <- "summary.ds"
  return(ans)
}
