# see predict.ds for documentation
#' @S3method predict rem
predict.rem <-
function(object,newdata=NULL,compute=FALSE,int.range=NULL,...)
#
# predict.rem - computes fitted values for p 
# 
# arguments:
#
# object     - rem model object
# newdata    - data for which predictions made (not currently used)
# compute    - if TRUE re-compute fitted values even if they are in model
# int.range  - integration range for variable range analysis
#
# return value:
#      list with 1 element:
#         fitted: vector of fitted detection probabilities - integral of p(y)/W
#
# Functions Used: predict.rem.fi and predict.ds 
#
#
{
   model=object
   if(!is.null(newdata))
   {
      compute=TRUE
      xmat=newdata
   }
   else
   xmat=model$mr$data
   xmat$distance=0
   ddfobj=model$ds$ds$aux$ddfobj
   if(ddfobj$type=="gamma")
   {
	   key.scale <- scalevalue(ddfobj$scale$parameters,ddfobj$scale$dm)
	   key.shape <- scalevalue(ddfobj$shape$parameters,ddfobj$shape$dm)
	   xmat$distance=rep(apex.gamma(key.scale,key.shape),2)
   }
   xmat$offsetvalue=0
   p.0=predict(model$mr,newdata=xmat,integrate=FALSE,compute=compute)$fitted
   if(is.null(newdata))
      pdot=predict(model$ds,esw=FALSE,compute=compute,int.range=int.range)$fitted
   else
      pdot=predict(model$ds,newdata=newdata[newdata$observer==1,],esw=FALSE,compute=compute,int.range=int.range)$fitted
   fitted=p.0*pdot
   if(is.null(newdata))
      names(fitted)=model$mr$mr$data$object[model$mr$mr$data$observer==1]
   else
      names(fitted)=newdata$object[newdata$observer==1]
   return(list(fitted=fitted))
}


