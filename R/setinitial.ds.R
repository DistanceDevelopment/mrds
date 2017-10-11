#' Set initial values for detection function based on distance sampling
#'
#' For a given detection function, it computes the initial values for the
#' parameters including scale and shape parameters and adjustment function
#' parameters if any.  If there are user-defined initial values only the
#' parameters not specified by the user are computed.
#'
#' @usage setinitial.ds(ddfobj, width, initial, point, left)
#'        sethazard(ddfobj, dmat, width, left)
#' @aliases setinitial.ds sethazard
#' @param ddfobj distance detection function object
#' @param width half-width of transect or radius of point count
#' @param left left truncation
#' @param initial list of user-defined initial values with possible elements
#'   scale,shape,adjustment
#' @param point if TRUE, point count data; otherwise, line transect data
#' @param dmat xmat from ddfobj
#' @return \item{scale}{vector of initial scale parameter values}
#'   \item{shape}{vector of initial shape parameter values}
#'   \item{adjustment}{vector of initial adjustment function parameter values}
#' @author Jeff Laake, David L Miller
#' @importFrom stats lm setNames
setinitial.ds <- function(ddfobj, width, initial, point, left){

  # storage
  if(ddfobj$type == "unif"){
    initialvalues <- list(scale=NULL, shape=NULL)
  }else{
    initialvalues <- list()
  }

  dmat <- ddfobj$xmat
  if(point){
    dmat$distance <- sqrt(dmat$distance)
  }
  point <- FALSE

  if(!any(is.na(initial)) & !is.list(initial)){
    stop("\ninitial values must be specified as a list with possible elements scale, shape, adjustments")
  }

  # Set parameters for hazard-rate
  if(ddfobj$type == "hr"){

    if(!all(is.na(initial)) && !is.null(initial$shape)){
      # use user-supplied value
      if(ncol(ddfobj$shape$dm) == length(initial$shape)){
        initialvalues$shape <- initial$shape
      }else{
        stop("Length of initial values for shape incorrect")
      }
    }else{
      # else find initial value
      initialvalues <- sethazard(ddfobj, dmat, width, left)
      if(ncol(ddfobj$shape$dm)>1){
        initialvalues$shape <- c(initialvalues$shape,
                                 rep(0, ncol(ddfobj$shape$dm)-1))
      }
    }

    # scale parameter for the hr
    if(!all(is.na(initial)) && !is.null(initial$scale)){
      if(ncol(ddfobj$scale$dm) == length(initial$scale)){
        initialvalues$scale <- initial$scale
      }else{
        stop("Length of initial values for scale incorrect")
      }
    }else{
        initialvalues$scale <- lm(eval(parse(text=paste(
                              "log(distance+width/1000)", ddfobj$scale$formula))),
                              data=dmat[dmat$detected==1, ])$coeff
    }
  }else{
    # otherwise...
    if(ddfobj$type!="unif"){
      if(!all(is.na(initial)) && !is.null(initial$scale)){
        if(ncol(ddfobj$scale$dm) == length(initial$scale)){
          initialvalues$scale <- initial$scale
        }else{
          stop("Length of initial values for scale incorrect")
        }
      }else{
        # Set scale parameters using Ramsey's approach of linear model
        # with log(distance)
        initialvalues <- list(scale=lm(eval(parse(text=paste(
                              "log(distance+width/1000)", ddfobj$scale$formula))),
                              data=dmat[dmat$detected==1, ])$coeff)

      }
    }

    # Set shape parameter values in a very cheesey way...
    if(!is.null(ddfobj$shape)){
      initialvalues$shape <- c(log(2), rep(0, ncol(ddfobj$shape$dm)-1))
    }
  }

  # Set initial values for the adjustment term parameters
  if(!is.null(ddfobj$adjustment)){
    initialvalues$adjustment <- rep(0, length(ddfobj$adjustment$order))
  }
  if(!all(is.na(initial)) && !is.null(initial$adjustment)){
    if(length(initialvalues$adjustment) == length(initial$adjustment)){
      initialvalues$adjustment <- initial$adjustment
    }else{
      stop("Length of initial values for adjustments incorrect")
    }
  }

  # okay, done?
  return(initialvalues)
}
