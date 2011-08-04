#' Predictions from distance sampling rem models
#' 
#' Predict detection probabilities values from a fitted
#' distance sampling rem model using either the original data or a new dataframe.
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
#' @aliases predict.rem 
#' @param object \code{ddf} model object
#' @param newdata new dataframe for prediction
#' @param compute if TRUE compute values and don't use the fitted values stored
#'   in the model object
#' @param int.range integration range for variable range analysis
#' @param \dots unspecified and unused arguments for S3 consistency
#' @export
#' @method predict rem
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
#' @seealso \code{\link{ddf}}, \code{\link{summary.rem}}
#' @keywords utility
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


