# see predict.ds for documentation
#' @S3method predict rem.fi
predict.rem.fi <-
function(object,newdata=NULL,compute=FALSE, int.range=NULL,integrate=FALSE,...)
#
# predict.rem.fi - computes fitted values for p 
#                    if model object already contains fitted values it uses those, otherwise it calculates them
# 
# arguments:
#
# object     - rem.fi model object
# newdata    - data for which predictions made 
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
       newdata=newdata[newdata$object %in% as.numeric(model$data$object),]
    }
    newdata$offsetvalue=0
	newdata2=newdata[newdata$observer==2,]
	newdata1=newdata[newdata$observer==1,]

    if(!integrate)
    {
	   p1 <- predict(model$mr,newdata1,type="response")
	   p2 <- predict(model$mr,newdata2,type="response")
	   fitted=p1+p2-p1*p2
	   if(is.null(newdata))
		   names(fitted)=model$mr$data$object[model$mr$data$observer==1]
	   else
		   names(fitted)=newdata$object[newdata$observer==1]
	   return(list(fitted=fitted,p1=p1,p2=p2))
   } else
   {
       left=model$meta.data$left
       formula=paste("~",as.character(model$mr$formula)[3],collapse="")
       if (class(model$mr)[1]=="gam")
         integral.numeric<-TRUE
       else
         integral.numeric <- is.linear.logistic(newdata,formula,length(coef(model$mr)),width)
       models<- list(g0model=formula,scalemodel=NULL,fullscalemodel=NULL)
       if(is.null(int.range)) 
          pdot.list <- pdot.dsr.integrate.logistic(width, width, coef(model$mr), newdata,integral.numeric, FALSE, models, rem=TRUE, point=point)
       else
          pdot.list <- pdot.dsr.integrate.logistic(int.range, width, coef(model$mr), newdata,integral.numeric, FALSE, models, rem=TRUE, point=point)
       if(left !=0) pdot.list$pdot <- pdot.list$pdot -
                       pdot.dsr.integrate.logistic(left, width, coef(model$mr), newdata,integral.numeric, FALSE, models, rem=TRUE, point=point)$pdot
       fitted <- pdot.list$pdot
   }
   if(is.null(newdata))
      names(fitted)=model$mr$data$object[model$mr$data$observer==1]
   else
      names(fitted)=newdata$object[newdata$observer==1]
   return(list(fitted=fitted))
}
