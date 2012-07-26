#' Predictions from distance sampling ds models
#' 
#' Predict detection probabilities (or esw) values from a fitted
#' distance sampling model using either the original data or a new dataframe.
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
#' @aliases predict.ds predict.ddf 
#' @param object \code{ddf} model object
#' @param newdata new dataframe for prediction
#' @param compute if TRUE compute values and don't use the fitted values stored
#'   in the model object
#' @param int.range integration range for variable range analysis
#' @param esw if TRUE, returns effective strip half-width (or effective
#'   detection radius for points) integral 0-W p(y)dy; otherwise it returns
#'   integral 0-W py)*pi(y) where pi(y)=1/W for lines and pi(y)=2r/W^2 for
#'   points.
#' @param \dots unspecified and unused arguments for S3 consistency
#' @export
#' @method predict ds
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
#' @seealso \code{\link{ddf}}, \code{\link{summary.ds}},
#'   \code{\link{plot.ds}}
#' @keywords utility
# Uses: integratedetfct, integratedetfct.logistic
predict.ds <- function(object,newdata=NULL,compute=FALSE,int.range=NULL,
                       esw=FALSE,...){
	model <- object
	ltmodel <- model$ds
	x <- ltmodel$aux$ddfobj$xmat   
	point <- ltmodel$aux$point
	width <- ltmodel$aux$width

  # If there are no fitted values present or compute is set TRUE or
  # a newdata frame has been used, then the predicted values must be computed.
  if(is.null(model$fitted) | compute | !is.null(newdata)){

    # Get model and extract the parameters. Note that in computing derivatives
    # for variances, the model parameters are perturbed and may not be at
    # mle values after fitting.
    fpar <- model$par
	  ddfobj <- ltmodel$aux$ddfobj
	  ddfobj <- assign.par(ddfobj,fpar)
	  doeachint <- ltmodel$aux$doeachint

    # Get integration ranges either from specified argument or from
    # values stored in the model.
    if(is.null(int.range)){
      if(is.null(ltmodel$aux$int.range)){
        int.range <- cbind(rep(0,nrow(ddfobj$xmat)),
                           rep(width,nrow(ddfobj$xmat)))
      }else{
        int.range <- ltmodel$aux$int.range
	      if(is.vector(int.range)){
          int.range <- cbind(rep(int.range[1],nrow(ddfobj$xmat)),
                             rep(int.range[2],nrow(ddfobj$xmat)))
        }
      }
    }

    # Extract other values from model object
    if(!is.null(newdata)){
	    if(!is.null(ddfobj$scale)){
		    zdim <- ncol(ddfobj$scale$dm)
		    znames <- colnames(ddfobj$scale$dm)
		    ddfobj$scale$dm <- setcov(newdata, as.formula(ddfobj$scale$formula))$cov
		    if(zdim != ncol(ddfobj$scale$dm) | 
           !all(znames==colnames(ddfobj$scale$dm)) ){
			    stop("fields or factor levels in newdata do not match data used in estimation model for scale model\n")
        }
      }
	  
      if(!is.null(ddfobj$shape)){
		    zdim <- ncol(ddfobj$shape$dm)
		    znames <- colnames(ddfobj$shape$dm)
		    ddfobj$shape$dm <- setcov(newdata, as.formula(ddfobj$shape$formula))$cov
		    if(zdim != ncol(ddfobj$shape$dm) | 
           !all(znames==colnames(ddfobj$shape$dm))){
			    stop("fields or factor levels in newdata do not match data used in estimation model for shape model\n")
        }
	    }
      # update xmat too
      datalist <- process.data(newdata,object$meta.data,check=FALSE)
      ddfobj$xmat <- datalist$xmat
    }

    # Compute integral of fitted detection function using either logistic or
    # non-logistic detection function.  Note that "logistic" is not currently
    # allowed as it has not been fully tested.

    # if(ftype=="logistic")
    #   int1=integratedetfct.logistic(x,ltmodel$model$scalemodel,width,
    #                         int.range,theta1,ltmodel$aux$integral.numeric,z)
    # else
  	ddfobj$cgftab <- tablecgf(ddfobj,width=width,standardize=TRUE, point=point)
    int1 <- integratepdf(ddfobj,select=rep(TRUE,nrow(ddfobj$xmat)),width=width,
		                     int.range=int.range,doeachint=doeachint,
                         standardize=TRUE,point=point)
    # int1=integratedetfct(ddfobj,select=rep(TRUE,nrow(ddfobj$xmat)),
    #           width=width,int.range=int.range,doeachint=doeachint,point=point)

    # If the predicted values don't need to be computed, then use the values 
    # in the model object (model$fitted) and change to integral (esw) values.
    # Note this needs to be checked to see if it works with variable ranges.

  }else{
    int1 <- model$fitted
  }

  # Compute either esw (int1) or p and store in fitted.
  if(esw){
	  if(!point){
      fitted <- int1*width
    }else{
	    fitted <- int1*pi*width^2
    }
  }else{
	   fitted <- int1	
  }

  # If there are no covariates and there is only one prediction, expand to 
  # a vector based on length of data object.
  if(length(fitted)==1){
    if(is.null(newdata)){
      fitted <- rep(fitted,length(x$object))
    }else{
      fitted <- rep(fitted,nrow(newdata))
    }
  }

  # If not a new dataframe, then assign names from data stored in the model 
  # object otherwise, use those from newdata.  Then return vector of values 
  # to calling frame.
  if(is.null(newdata)){
    names(fitted) <- x$object
  }else{
    names(fitted) <- newdata$object
  }
  
  return(list(fitted=fitted))
}
