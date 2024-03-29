#' Summary of distance detection function model object
#'
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error.
#'
#' The argument \code{N} is used to suppress computation of
#' abundance and average detection probability in calls to summarize the
#' \code{ds} and either \code{io.fi} or \code{trial.fi} for summaries of
#' \code{io} and \code{trial} objects respectively which are composed of a
#' \code{ds} model object and a mark-recapture model object. The corresponding
#' print function is called to print the summary results.
#'
#' @export
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param N if TRUE, computes abundance in covered (sampled) region
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.ds <- function(object, se=TRUE, N=TRUE, ...){

  model <- object
  avgp <- function(model, pdot, ...){return(pdot)}

  ddfobj <- model$ds$aux$ddfobj

  ans <- list()

  # was monotonicity enforced?
  ans$mono <- model$ds$aux$mono
  # strict monotonicity?
  ans$mono.strict <- model$ds$aux$mono.strict

  # Number of observations
  ans$n <- length(ddfobj$xmat$distance)

  # Set the key function type
  ans$key <- ddfobj$type

  # Parameter estimates and se for detection function
  # se is included as part of the objects, see coef.ds
  # for details
  coeff <- coef(model)
  if(!is.null(coeff)){
    # Scale Coefficient
    ans$coeff$key.scale <- coeff$scale

    # Hazard shape parameter
    if(ans$key%in%c("gamma", "hr", "th1", "th2", "tpn")){
     ans$coeff$key.shape <- coeff$exponent
    }

    # Adjustment term parameter(s)
    # See coef.ds() on how this is returned
    # This is a vector remember, so if you are using this
    # you need to take that into account.
    if(!is.null(coeff$adjustment)){
      ans$adjustment <- ddfobj$adjustment
      ans$coeff$adj.order <- model$adj.order
      ans$coeff$adj.parm <- coeff$adjustment
    }
  }else{
    ans$coeff <- NULL
  }

  # AIC
  ans$aic <- model$criterion
  
  # Optimisation
  ans$optimise <- model$optimise

  # Truncation distances, left and right
  ans$width <- model$meta.data$width
  ans$left <- model$meta.data$left

  ans$average.p <- ans$n/model$Nhat

  # find se of average p and Nhat (in covered area)
  if(se &!is.null(ans$coeff)){
    se.obj <- calc.se.Np(model, avgp, ans$n, ans$average.p)
  }

  if(N){
    ans$Nhat <- model$Nhat
    if(se&!is.null(ans$coeff)){
      ans$Nhat.se <- se.obj$Nhat.se
    }
  }

  if(se & !is.null(ans$coeff)){
    ans$average.p.se <- se.obj$average.p.se
  }

  # save transect type
  ans$transect <- c("line", "point")[model$meta.data$point+1]

  # flag if the integration ranges were set
  ans$int.range <- FALSE
  if(is.matrix(object$meta.data$int.range) &&
     nrow(unique(object$meta.data$int.range)) > 1){
    ans$int.range <- TRUE
  }

  class(ans) <- "summary.ds"
  return(ans)
}
