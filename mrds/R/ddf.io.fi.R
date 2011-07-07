

#' Mark-Recapture Analysis of Independent Observer Configuration with Full
#' Independence
#' The mark-recapture data derived from an independent observer distance
#' sampling survey can be used to derive conditional detection functions
#' (p_j(y)) for both observers (j=1,2).  They are conditional detection
#' functions because detection probability for observer j is based on seeing or
#' not seeing observations made by observer 3-j. Thus, p_1(y) is estimated by
#' p_1|2(y).  If detections by the observers are independent (full
#' independence) then p_1(y)=p_1|2(y),p_2(y)=p_2|1(y) and for the union, full
#' independence means that p(y)=p_1(y) + p_2(y) - p_1(y)*p_2(y) for each
#' distance y.  In fitting the detection functions the likelihood given by eq
#' 6.8 and 6.16 in Laake and Borchers (2004) is used. That analysis does not
#' require the usual distance sampling assumption that perpendicular distances
#' are uniformly distributed based on line placement that is random relative to
#' animal distribution.  However, that assumption is used in computing
#' predicted detection probability which is averaged based on a uniform
#' distribution (see eq 6.11 of Laake and Borchers 2004).
#' 
#' For a complete description of each of the calling arguments, see
#' \code{\link{ddf}}.  The argument \code{model} in this function is the same
#' as \code{mrmodel} in \code{ddf}.  The argument \code{dataname} is the name
#' of the dataframe specified by the argument \code{data} in \code{ddf}. The
#' arguments \code{control},\code{meta.data},and \code{method} are defined the
#' same as in \code{ddf}.
#' 
#' @method ddf io.fi
#' @param model mark-recapture model specification
#' @param data analysis dataframe
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call original function call used to call \code{ddf}
#' @param method analysis method; only needed if this function called from
#'   \code{ddf.io}
#' @return result: an io.fi model object
#' @export
#' @author Jeff Laake
#' @seealso
#'   \code{\link{ddf.io}},\code{\link{summary.io.fi}},\code{\link{coef.io.fi}},\code{\link{plot.io.fi}},
#'   \code{\link{gof.io.fi}},\code{\link{io.glm}}
#' @references Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#' @keywords Statistical Models
ddf.io.fi <-
function(model,data,meta.data=list(),control=list(),call="",method)
{
# 
# ddf.io.fi
#
# Fits double-observer data with io configuration and full independence (L_omega only).
#
# Arguments:
#
#  model     - mr model object
#  data      - dataframe
#  meta.data - list containing settings controlling data structure
#  control   - list containing settings controlling model fitting
#  call      - call used for ddf
#  method    - io or io.fi (if io this is being called as part of fitting) 
#
#  Functions used: assign.default.values, process.data, create.model.frame
#                  ioglm, predict(predict.io.fi), NCovered (NCovered.io.fi)
#
#  Return value: model object of class "io.fi"
#
#  NOTE: gams are only partially implemented
# 
# The following are dummy glm and gam functions that are defined here to provide the
# list of arguments for use in the real glm/gam functions.  These dummy functions are
# removed after they are used so the real ones can be used in the model fitting.
#
glm=function(formula,link="logit")
{
if(class(formula)!="formula")
{
   if(class(try(as.formula(formula)))!="formula")
      stop("Invalid formula")
}
else
   formula=paste(as.character(formula),collapse="") 
if(class(link)=="function") link=substitute(link)
link=match.arg(link,c("logit"))
return(list(fct="glm",formula=formula,link=substitute(link)))
}
gam=function(formula,link="logit")
{
if(class(formula)!="formula")
{
   if(class(try(as.formula(formula)))!="formula")
      stop("Invalid formula")
}
else
   formula=paste(as.character(formula),collapse="") 
if(class(link)=="function") link=substitute(link)
link=match.arg(link,c("logit"))
return(list(fct="gam",formula=formula,link=substitute(link)))
}
#
#  Save current user options and then set design contrasts to treatment style
#
   save.options<-options()
   options(contrasts=c("contr.treatment","contr.poly"))
#
# Set up meta data values
#
	meta.data=assign.default.values(meta.data, left=0, width=NA, binned=FALSE, int.range=NA, strict=TRUE,
			nonmono=FALSE,fdebug=0,engine="optim",point=FALSE)
#
# Set up control values
#
   control=assign.default.values(control,showit = FALSE, doeachint=FALSE, estimate=TRUE,refit=TRUE,nrefits=25,
                                       initial = NA, lowerbounds = NA, upperbounds = NA)
#
# Assign model values; this uses temporarily defined functions glm and gam
#
   modpaste=paste(model)
   modelvalues=try(eval(parse(text=modpaste[2:length(modpaste)])))
   if(class(modelvalues)=="try-error")stop("Invalid model specification: ",model)
   rm(glm,gam)
#
#  Process data if needed  
#
   if(is.data.frame(data))
   {
      data.list=process.data(data,meta.data)
      meta.data=data.list$meta.data
      xmat=data.list$xmat
   }
   else
   {                  
      xmat=data
   }
#
#  Setup default breaks
#
   if(meta.data$binned)
      meta.data$breaks=c(max(0,min(as.numeric(levels(as.factor(xmat$distbegin))))),as.numeric(levels(as.factor(xmat$distend))))
#
#  Create result list with some arguments
#
   result=list(call=call,data=data,model=model,meta.data=meta.data,
                    control=control,method="io.fi")
   class(result)=c("io.fi","ddf")
#
#  Create formula and model frame (if not GAM)
#
   xmat$offsetvalue <- rep(0,dim(xmat)[1])
   model.formula=paste("detected",modelvalues$formula)
   GAM=FALSE
   if(modelvalues$fct=="gam") 
   {
      GAM=TRUE
   } else
      xmat=create.model.frame(xmat,as.formula(model.formula),meta.data)
   model.formula=as.formula(paste(model.formula,"+offset(offsetvalue)"))
#
#  Fit the conditional detection functions using io.glm 
#
   result$mr <- io.glm (xmat,model.formula,GAM=GAM)
   if(GAM)result$mr$data=xmat
#
#  Compute the L_omega portion of the likelihood value, AIC and hessian
#
   cond.det=predict(result)
   result$par <- coef(result$mr)  
   npar=length(result$par)
   p1=cond.det$p1
   p2=cond.det$p2
   p.c.omega=p1^result$mr$data$detected[result$mr$data$observer==1]*
		   (1-p1)^(1-result$mr$data$detected[result$mr$data$observer==1])*
		   p2^result$mr$data$detected[result$mr$data$observer==2]*
		   (1-p2)^(1-result$mr$data$detected[result$mr$data$observer==2]) 
   result$lnl <- sum(log(p.c.omega)) - sum(log(cond.det$fitted))      
   if(GAM) 
      result$hessian <- result$mr$Vp
   else
      result$hessian <- solve(summary(result$mr)$cov.unscaled) 

#
#  If this is method=io.fi then compute lnl for L_y and add to L_omega before
#  computing AIC. Note the code in predict is currently specific to logit link.    
#
   if(method=="io.fi")
   {
#
#     Compute fitted values - integral of p(x)/width    
#
	  result$fitted=predict(result,integrate=TRUE)$fitted
      result$Nhat=NCovered(result$par,result)
	  distances=result$data$distance[result$data$observer==1]
	  if(meta.data$point)
		 result$lnl<- result$lnl + sum(log(cond.det$fitted*2*distances/meta.data$width^2)) - sum(log(result$fitted))
	  else
         result$lnl<- result$lnl + sum(log(cond.det$fitted/meta.data$width)) - sum(log(result$fitted))
   }
   result$criterion<- -2*result$lnl + 2*npar
#
# Restore user options
#
   options(save.options)
#
# Return result
#
   return(result)
}
