#' Predictions from distance sampling trial.fi models
#' 
#' Predict detection probabilities values from a fitted
#' distance sampling trial.fi model using either the original data or a new dataframe.
#' 
#' The first 4 arguments are the same in each predict function.  The latter 2
#' are specific to certain functions. The effective strip half-width (esw) is
#' the integral of the fitted detection function over the range of the sampled
#' area (either 0 to W or the specified \code{int.range}).  The predicted
#' detection probability is the average probability which is simply the
#' integral divided by the distance range.  The fitted detection probabilities
#' are stored in the \code{model} object and these are used unless
#' \code{compute=TRUE} or \code{newdata} is specified. \code{compute=TRUE} is
#' used to estimate numerical derivatives for use in delta method
#' approximations to the variance.  For \code{method="io.fi" or ="trial.fi"} if
#' \code{integrate=FALSE}, \code{predict} returns the value of the conditional
#' detection probability and if \code{integrate=TRUE}, it returns the average
#' conditional detection probability by integrating over x(distance) with
#' respect to a uniform distribution.
#' 
#' @aliases predict.trial.fi 
#' @param object \code{ddf} model object
#' @param newdata new dataframe for prediction
#' @param compute if TRUE compute values and don't use the fitted values stored
#'   in the model object
#' @param int.range integration range for variable range analysis
#' @param integrate if TRUE compute integral over p(y)*pi(y); otherwise just
#'   compute p(y)
#' @param \dots unspecified and unused arguments for S3 consistency
#' @export
#' @method predict trial.fi
#' @return For all but the exceptions below,the value is a list with a single
#'   element: \tabular{ll}{ \code{fitted} \tab vector of average detection
#'   probabilities or esw values for each observation in the original data or
#'   \code{newdata} \cr }
#' 
#' For \code{predict.io.fi},\code{predict.trial.fi},\code{predict.rem.fi} with
#'   integrate=TRUE, he value is a list with the elements: \tabular{ll}{
#'   \code{fitted} \tab vector of integrated (average) detection probabilities
#'   for each observation in the original data or \code{newdata} \cr }
#' 
#' For \code{predict.io.fi}, \code{predict.trial.fi}, or \code{predict.rem.fi}
#'   with \code{integrate=FALSE}, the value is a list with the following
#'   elements: \tabular{ll}{ \code{fitted} \tab p(y) values \cr \code{p1} \tab
#'   p_1|2(y) (conditional detection probability for observer 1) \cr \code{p2}
#'   \tab p_2|1(y) (conditional detection probability for observer 2) \cr
#'   \code{fitted} \tab p_.(y)=p_1|2(y)+p_2|1(y)-p_1|2(y)*p_2|1(y) (conditional
#'   detection probability of being seen by either observer) \cr }
#' @note Each function is called by the generic function \code{predict} for the
#'   appropriate \code{ddf} model object.  They can be called directly by the
#'   user, but it is typically safest to use \code{predict} which calls the
#'   appropriate function based on the type of model.
#' @author Jeff Laake
#' @seealso \code{\link{ddf}}, \code{\link{summary.trial.fi}},
#'   \code{\link{plot.trial.fi}}
#' @keywords utility
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
