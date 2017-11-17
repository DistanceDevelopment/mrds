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

  ftype <- ddfobj$type
  if(ftype == "unif"){
    initialvalues <- list(scale=NULL, shape=NULL)
  }

  dmat <- ddfobj$xmat
  # unclear of the origin of this code but somehow it messes everything up
  #if(point){
  #  dmat$distance <- sqrt(dmat$distance)
  #}
  point <- FALSE

  # Set shape parameters for special case of cds hazard function
  if(ftype == "hr"){
    initialvalues <- sethazard(ddfobj, dmat, width, left)
    if(ncol(ddfobj$shape$dm)>1){
      initialvalues$shape <- c(initialvalues$shape,
                               rep(0, ncol(ddfobj$shape$dm)-1))
    }
    if(ncol(ddfobj$scale$dm)>1){
      initialvalues$scale <- lm(eval(parse(text=paste(
                        "log(distance+width/1000)", ddfobj$scale$formula))),
                            data=dmat[dmat$detected==1, ])$coeff
    }
  }else{
    # Set scale parameters using Ramsey's approach of linear model
    # with log(distance)
    if(ftype!="unif"){
      initialvalues <- list(scale=lm(eval(parse(text=paste(
                            "log(distance+width/1000)", ddfobj$scale$formula))),
                            data=dmat[dmat$detected==1, ])$coeff)
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

  if(!any(is.na(initial))){
    if(!is.list(initial)){
      stop("\ninitial values must be specified as a list with possible elements scale, shape, adjustments")
    }
    if(!is.null(initial$shape)){
      if(length(initialvalues$shape) == length(initial$shape)){
        initialvalues$shape <- setNames(initial$shape,
                                        names(initialvalues$shape))
      }else{
        stop("Length of initial values for shape incorrect")
      }
    }
    if(!is.null(initial$scale)){
      if(length(initialvalues$scale) == length(initial$scale)){
        initialvalues$scale <- setNames(initial$scale,
                                        names(initialvalues$scale))
      }else{
        stop("Length of initial values for scale incorrect")
      }
    }

    if(!is.null(initial$adjustment)){
      if(length(initialvalues$adjustment) == length(initial$adjustment)){
        initialvalues$adjustment <- setNames(initial$adjustment,
                                             names(initialvalues$adjustment))
      }else{
        stop("Length of initial values for adjustments incorrect")
      }
    }
  }
  return(initialvalues)
}
