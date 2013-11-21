# see predict.ds for documentation
#' @S3method predict trial
predict.trial <- function(object,newdata=NULL,compute=FALSE,int.range=NULL,...){
#
# predict.trial - computes fitted values for p 
# 
# arguments:
#
# object     - trial model object
# newdata    - data for which predictions made (only used to change distance=0)
# compute    - if TRUE re-compute fitted values even if they are in model
# int.range  - integration range for variable range analysis
# new        - if TRUE, newdata is truly new data and not just a subset for trials
#
# return value
#  list with 1 element - vector of fitted detection probabilities 
# 
# Functions Used: predict.trial.fi and predict.ds 
#    
   model=object
   if(is.null(newdata))
   {
      xmat=model$data
      xmat=xmat[xmat$observer==1&xmat$object%in%as.numeric(names(model$fitted)),]
   }
   else
   {
      xmat=newdata
      compute=TRUE
   }
   if(!compute)
      return(list(fitted=model$fitted))
   xmat$distance=0
   ddfobj <- model$ds$ds$aux$ddfobj
   if(ddfobj$type=="gamma")
   {
	   key.scale <- scalevalue(ddfobj$scale$parameters,ddfobj$scale$dm)
	   key.shape <- scalevalue(ddfobj$shape$parameters,ddfobj$shape$dm)
	   xmat$distance=as.vector(apex.gamma(key.scale,key.shape))
   }    
   if(is.null(newdata))
      pdot=predict(model$ds,esw=FALSE,compute=compute,int.range=int.range)$fitted
   else
      pdot=predict(model$ds,newdata=newdata,esw=FALSE,compute=compute,int.range=int.range)$fitted
   fitted=predict(model$mr$mr,newdata=xmat,type="response")*pdot
   if(is.null(newdata))
      names(fitted)=names(model$fitted)
   else
      names(fitted)=newdata$object[newdata$observer==1]
   return(list(fitted=fitted))
}
