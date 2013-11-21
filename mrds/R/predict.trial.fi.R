#  see predict.ds for documentation
#' @S3method predict trial.fi
predict.trial.fi <-
function(object,newdata=NULL,compute=FALSE, int.range=NULL,integrate=FALSE,...)
#
# predict.trial.fi - computes fitted values for p 
#                    if model object already contains fitted values it uses those, otherwise it calculates them
# 
# arguments:
#
# object     - trial.fi model object
# newdata    - data for which predictions made - only used to change value of distance to 0
# compute    - if TRUE re-compute fitted values even if they are in model
#              (not currently used)
# int.range  - integration range for variable range analysis
# integrate  - if TRUE compute integral of p(x) otherwise just p(x)
# new        - if TRUE, newdata is truly new data and not just a subset for trials
#
# return value
#  list with 1 element - vector of fitted detection probabilities 
# 
# Functions Used: pdot.dsr.integrate.logistic, is.linear.logistic, predict.glm (could also use predict.gam eventually) 
{
   model=object
   point=model$meta.data$point
   width=model$meta.data$width
   if(is.null(newdata))
   {
       newdata=model$data
       newdata=newdata[newdata$object %in% as.numeric(names(model$fitted)),]
   }
   if(!integrate)
   {
		p1<-predict(model$mr,newdata,type="response")
		p2<-0
		fitted=p1
		if(is.null(newdata))
			names(fitted)=model$mr$data$object[model$mr$data$observer==1]
		else
			names(fitted)=newdata$object[newdata$observer==1]
		return(list(p1=p1,p2=0,fitted=fitted))
	}
   else
   {
       left=model$meta.data$left
       formula=paste("~",as.character(model$mr$formula)[3],collapse="")
       if (class(model$mr)[1]=="gam")
         integral.numeric<-TRUE
       else
         integral.numeric <- is.linear.logistic(newdata,formula,length(coef(model$mr)),width)
       models<- list(g0model=formula,scalemodel=NULL,fullscalemodel=NULL)
       if(is.null(int.range)) 
          pdot.list <- pdot.dsr.integrate.logistic(width, width, coef(model$mr), newdata,integral.numeric, TRUE, models, point=point)
       else
          pdot.list <- pdot.dsr.integrate.logistic(int.range, width, coef(model$mr), newdata,integral.numeric, TRUE, models, point=point)
       if(left !=0) pdot.list$pdot <- pdot.list$pdot -
                       pdot.dsr.integrate.logistic(left, width, coef(model$mr), newdata,integral.numeric, TRUE, models, point=point)$pdot
	   fitted=pdot.list$pdot
   }
   if(is.null(newdata))
      names(fitted)=model$mr$data$object[model$mr$data$observer==1]
   else
      names(fitted)=newdata$object[newdata$observer==1]
   return(list(fitted=fitted))
}
