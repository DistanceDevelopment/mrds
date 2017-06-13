#' Predictions from \code{mrds} models
#'
#' Predict detection probabilities (or effective strip widths/effective areas of detection) from a fitted distance sampling model using either the original data (i.e. "fitted" values) or using new data.
#'
#' The first 4 arguments are the same in each predict function.  The latter 2 are specific to certain functions. For line transects, the effective strip half-width (\code{esw=TRUE}) is the integral of the fitted detection function over either 0 to W or the specified \code{int.range}.  The predicted detection probability is the average probability which is simply the integral divided by the distance range. For point transect models, \code{esw=TRUE} calculates the effective area of detection (commonly referred to as "nu", this is the integral of \code{2/width^2 * rg(r)}.
#'
#' Fitted detection probabilities are stored in the \code{model} object and these are returned unless \code{compute=TRUE} or \code{newdata} is specified. \code{compute=TRUE} is used to estimate numerical derivatives for use in delta method approximations to the variance.
#'
#' For \code{method="io.fi"} or \code{method="trial.fi"} if \code{integrate=FALSE}, \code{predict} returns the value of the conditional detection probability and if \code{integrate=TRUE}, it returns the average conditional detection probability by integrating over x (distance) with respect to a uniform distribution.
#'
#' Note that the ordering of the returned results when no new data is supplied (the "fitted" values) will not necessarily be the same as the data supplied to \code{\link{ddf}}, the data (and hence results from \code{predict}) will be sorted by object ID (\code{object}) then observer ID (\code{observer}).
#'
#' @aliases predict predict.ds predict.ddf predict.io predict.io.fi predict.trial predict.trial.fi predict.rem predict.rem.fi
#' @param object \code{ddf} model object.
#' @param newdata new \code{data.frame} for prediction, this must include a column called "\code{distance}".
#' @param compute if \code{TRUE} compute values and don't use the fitted values stored in the model object.
#' @param int.range integration range for variable range analysis; either vector or 2 column matrix.
#' @param esw if \code{TRUE}, returns effective strip half-width (or effective area of detection for point transect models) integral from 0 to the truncation distance (\code{width}) of \eqn{p(y)dy}; otherwise it returns the integral from 0 to truncation width of \eqn{p(y)\pi(y)} where \eqn{\pi(y)=1/w} for lines and \eqn{\pi(y)=2r/w^2} for points.
#' @param integrate for \code{*.fi} methods, see Details below.
#' @param \dots for S3 consistency
#' @usage \method{predict}{ds}(object, newdata, compute=FALSE, int.range=NULL, esw=FALSE, ...)
#'        \method{predict}{io.fi}(object, newdata, compute=FALSE, int.range=NULL, integrate=FALSE, ...)
#'        \method{predict}{io}(object, newdata, compute=FALSE, int.range=NULL, ...)
#'        \method{predict}{trial}(object, newdata, compute=FALSE, int.range=NULL, ...)
#'        \method{predict}{trial.fi}(object, newdata, compute=FALSE, int.range=NULL, integrate=FALSE, ...)
#'        \method{predict}{rem}(object, newdata, compute=FALSE, int.range=NULL, ...)
#'        \method{predict}{rem.fi}(object, newdata, compute=FALSE, int.range=NULL, integrate=FALSE, ...)
#' @return For all but the exceptions below, the value is a list with a single element: \code{fitted}, a vector of average detection probabilities or esw values for each observation in the original data or\code{newdata}
#'
#' For \code{predict.io.fi},\code{predict.trial.fi},\code{predict.rem.fi} with \code{integrate=TRUE}, the value is a list with one element: \code{fitted}, which is a vector of integrated (average) detection probabilities for each observation in the original data or \code{newdata}.
#'
#' For \code{predict.io.fi}, \code{predict.trial.fi}, or \code{predict.rem.fi} with \code{integrate=FALSE}, the value is a list with the following elements:
#'  \describe{
#'    \item{\code{fitted}}{\eqn{p(y)} values}
#'    \item{\code{p1}}{\eqn{p_{1|2}(y)}, conditional detection probability for observer 1}
#'    \item{\code{p2}}{\eqn{p_{2|1}(y)}, conditional detection probability for observer 2}
#'    \item{\code{fitted}}{\eqn{p_.(y)=p_{1|2}(y)+p_{2|1}(y)-p_{1|2}(y)*p_{2|1}(y)}, conditional detection probability of being seen by either observer}}
#'
#' @note Each function is called by the generic function \code{predict} for the appropriate \code{ddf} model object.  They can be called directly by the user, but it is typically safest to use \code{predict} which calls the appropriate function based on the type of model.
#' @author Jeff Laake, David L Miller
#' @seealso \code{\link{ddf}}, \code{\link{summary.ds}}, \code{\link{plot.ds}}
#' @keywords utility
#' @export
#' @importFrom stats predict
# Uses: integratedetfct, integratedetfct.logistic
predict.ds <- function(object, newdata=NULL, compute=FALSE, int.range=NULL,
                       esw=FALSE, ...){
  model <- object
  ltmodel <- model$ds
  x <- ltmodel$aux$ddfobj$xmat
  point <- ltmodel$aux$point
  width <- ltmodel$aux$width
  left <- model$meta.data$left
  ddfobj <- ltmodel$aux$ddfobj

  # Get integration ranges either from specified argument or from
  # values stored in the model.
  if(is.null(int.range)){
    if(is.null(newdata)){
      nr <- nrow(ddfobj$xmat)
    }else{
      nr <- nrow(newdata)
    }

    if(is.null(ltmodel$aux$int.range)){
      int.range <- cbind(rep(0, nr), rep(width, nr))
    }else{
      int.range <- ltmodel$aux$int.range
      if(is.vector(int.range)){
        int.range <- cbind(rep(int.range[1], nr),
                           rep(int.range[2], nr))
      #}else if(nrow(int.range) == (nrow(x)+1)){
      #int.range <- int.range[2:nrow(int.range), , drop=FALSE]
      }
    }
  }

  # If there are no fitted values present or compute is set TRUE or
  # a newdata frame has been used, then the predicted values must be computed.
  if(is.null(model$fitted) | compute | !is.null(newdata)){

    # Get model and extract the parameters. Note that in computing derivatives
    # for variances, the model parameters are perturbed and may not be at
    # mle values after fitting.
    fpar <- model$par
    ddfobj <- assign.par(ddfobj, fpar)

    # Extract other values from model object
    if(!is.null(newdata)){

      newdata_save <- newdata

      # get the data in the model
      model_dat <- model$data

      # do this for both scale and shape parameters
      for(df_par in c("scale", "shape")){
        # if that parameter exists...
        if(!is.null(ddfobj[[df_par]])){
          # save the column names from the design matrix
          znames <- colnames(ddfobj[[df_par]]$dm)

          # pull out the columns in the formula and the distances column
          fvars <- all.vars(as.formula(model$ds$aux$ddfobj[[df_par]]$formula))

          if(!all(fvars %in% colnames(newdata))){
            stop("columns in `newdata` do not match those in fitted model\n")
          }

          model_dat <- model_dat[,c("distance", fvars), drop=FALSE]

          # setup the covariate matrix, using the model data to ensure that
          # the levels are right
          newdata <- rbind(model_dat,
                           newdata_save[, c("distance", fvars), drop=FALSE])
          dm <- setcov(newdata, as.formula(ddfobj[[df_par]]$formula))

          # now check that the column names are the same for the model
          # and prediction data matrices
          if(!identical(colnames(dm), znames) || any(is.na(newdata))){
            stop("fields or factor levels in `newdata` do not match data used in fitted model\n")
          }

          # get only the new rows for prediction
          dm <- dm[(nrow(model_dat)+1):nrow(dm),,drop=FALSE]
          # assign that!
          ddfobj[[df_par]]$dm <- dm

        }
      }

      # handle data setup for uniform key case
      if(ddfobj$type == "unif"){
        model_dat <- model_dat[, "distance", drop=FALSE]
        newdata <- rbind(model_dat,
                         newdata_save[, "distance", drop=FALSE])
        dm <- setcov(newdata, ~1)
        dm <- dm[(nrow(model_dat)+1):nrow(dm), , drop=FALSE]
      }

      # get the bins when you have binned data
      # use the breaks specified in the model!
      if(model$meta.data$binned){
        newdata <- create.bins(newdata, model$meta.data$breaks)
      }

      # update xmat too
      datalist <- process.data(newdata, object$meta.data, check=FALSE)
      ddfobj$xmat <- datalist$xmat[(nrow(model_dat)+1):nrow(datalist$xmat),,drop=FALSE]
      # reset newdata to be the right thing
      newdata <- newdata[(nrow(model_dat)+1):nrow(newdata), , drop=FALSE]

    }

    # Compute integral of fitted detection function using either logistic or
    # non-logistic detection function.  Note that "logistic" is not currently
    # allowed as it has not been fully tested.

    # if(ftype=="logistic")
    #   int1=integratedetfct.logistic(x,ltmodel$model$scalemodel,width,
    #                         int.range,theta1,ltmodel$aux$integral.numeric,z)
    # else
    int1 <- integratepdf(ddfobj, select=rep(TRUE, nrow(ddfobj$xmat)),
                         width=width, int.range=int.range, standardize=TRUE,
                         point=point, left=left, doeachint=TRUE)
  }else{
    # If the predicted values don't need to be computed, then use the values
    # in the model object (model$fitted) and change to integral (esw) values.
    # Note this needs to be checked to see if it works with variable ranges.
    int1 <- model$fitted
  }

  # Compute either esw (int1) or p and store in fitted.
  if(esw){
    if(!point){
      fitted <- int1*(int.range[,2]-int.range[,1])
    }else{
      fitted <- int1*pi*width^2
    }
  }else{
     fitted <- int1
  }

  # If there are no covariates and there is only one prediction, expand to
  # a vector based on length of data object.
  if(length(fitted)==1){
    if(is.null(newdata)){
      fitted <- rep(fitted, length(x$object))
    }else{
      fitted <- rep(fitted, nrow(newdata))
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
