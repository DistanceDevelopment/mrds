#' Predictions from distance sampling io.fi models
#' 
#' Predict detection probabilities values from a fitted
#' distance sampling io.fi model using either the original data or a new dataframe.
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
#' @aliases predict.io.fi 
#' @param object \code{ddf} model object
#' @param newdata new dataframe for prediction
#' @param compute if TRUE compute values and don't use the fitted values stored
#'   in the model object
#' @param int.range integration range for variable range analysis
#' @param integrate if TRUE compute integral over p(y)*pi(y); otherwise just
#'   compute p(y)
#' @param \dots unspecified and unused arguments for S3 consistency
#' @export
#' @method predict io.fi
#' @return For all but the exceptions below,the value is a list with a single
#'   element: \tabular{ll}{ \code{fitted} \tab vector of average detection
#'   probabilities or esw values for each observation in the original data or
#'   \code{newdata} \cr }
#' 
#' For \code{predict.io.fi},\code{predict.trial.fi},\code{predict.rem.fi} with
#'   integrate=TRUE, the value is a list with the elements: \tabular{ll}{
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
#' @seealso \code{\link{ddf}}, \code{\link{summary.io.fi}},
#'   \code{\link{plot.io.fi}}
#' @keywords utility
predict.io.fi <- function(object,newdata=NULL,compute=FALSE, int.range=NULL, 
                          integrate=FALSE,...){
# Functions Used: pdot.dsr.integrate.logistic, is.linear.logistic, 
#                 predict.glm (could also use predict.gam eventually) 
  model <- object
  width <- model$meta.data$width
  point <- model$meta.data$point
  if(is.null(newdata)){
    newdata <- model$mr$data
  }
  newdata$offsetvalue <- 0
  GAM <- FALSE
  if("gam" %in% class(model$mr)){
    GAM <- TRUE
  }
  if(!integrate){
    fitted <- predict(model$mr,newdata,type="response")
    p1 <- fitted[model$mr$data$observer==1]
    p2 <- fitted[model$mr$data$observer==2]
    fitted <- p1+p2-p1*p2
    if(is.null(newdata)){
       names(fitted) <- model$mr$data$object[model$mr$data$observer==1]
    }else{
       names(fitted)=newdata$object[newdata$observer==1]
    }
    return(list(fitted=fitted,p1=p1,p2=p2))
  }else{
    left <- model$meta.data$left
    formula <- paste("~",as.character(model$mr$formula)[3],collapse="")
    if("gam" %in% class(model$mr)){
      integral.numeric <- TRUE
    }else{
      integral.numeric <- is.linear.logistic(newdata,formula,
                                             length(coef(model$mr)),width)
    }
    models <- list(g0model=formula,scalemodel=NULL,fullscalemodel=NULL)

    # now int.range is a vector with lower and upper bounds
    if(is.null(int.range)){
      pdot.list <- pdot.dsr.integrate.logistic(width,width, model$mr$coef, 
                     newdata,integral.numeric, FALSE, models,GAM, point=point)
    }else{
      pdot.list <- pdot.dsr.integrate.logistic(int.range,width, model$mr$coef, 
                       newdata,integral.numeric, FALSE, models,GAM, point=point)
    }

    if(left !=0){
      pdot.list$pdot <- pdot.list$pdot -
                    pdot.dsr.integrate.logistic(left, width, model$mr$coef, 
                                   newdata,integral.numeric, FALSE, models,GAM,
                                   point=point)$pdot
    }

    fitted <- pdot.list$pdot
    if(is.null(newdata)){
       names(fitted) <- model$mr$data$object[model$mr$data$observer==1]
    }else{
       names(fitted) <- newdata$object[newdata$observer==1]
    }

    return(list(fitted=fitted))
  }
}
