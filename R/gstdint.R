#' Integral of pdf of distances
#'
#' Computes the integral of distpdf with scale=1 (stdint=TRUE) or specified
#' scale (stdint=FALSE)
#'
#' @param x lower,upper value for integration
#' @param ddfobj distance detection function specification
#' @param index specific data row index
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count(TRUE) or line
#'   transect(FALSE)
#' @param stdint if TRUE, scale=1 otherwise specified scale used
#' @return vector of integral values of detection function
#' @note This is an internal function that is not intended to be invoked
#'   directly.
#' @author Jeff Laake and David L Miller
#' @keywords utility
gstdint <- function (x, ddfobj, index=NULL,select=NULL,width,
                      standardize=TRUE, point=FALSE, stdint=TRUE){

  ## NB this is not the integral of the PDF OR the detection function but
  ##    rather the integral of:
  ##       g(x)/w for line transects
  ##       2*r*g(r)/width^2

  # when we have half-normal, key only use the exact analytic expression
  # for the integral using the error function/analytic expression
  if(ddfobj$type=="hn" & is.null(ddfobj$adjustment)){

    # Set of observations for computation of detection function can
    # be specified with logical (select) and numeric (index) values.
    # Either or both can be specified although the latter is unlikely.
    if(is.null(select)){
      # use all
      if(is.null(index)){
        scale.dm <- ddfobj$scale$dm
      }else{
        # use only those with specific indices
        scale.dm <- ddfobj$scale$dm[index,,drop=FALSE]
      }
    }else{
      # Use those with select=TRUE
      if(is.null(index)){
        scale.dm <- ddfobj$scale$dm[select,,drop=FALSE]
      }else{
        # use the numeric index within those with select=TRUE
        scale.dm <- ddfobj$scale$dm[select,,drop=FALSE][index,,drop=FALSE]
      }
    }

    key.scale <- scalevalue(ddfobj$scale$parameters,scale.dm)
    if(point){
      # analytic expression for integral of 2*r*g(r)/width^2 when
      #  g(r) is half-normal
      int <- (2*(key.scale^2*exp(-x[1]^2/(2*key.scale^2))-
                 key.scale^2*exp(-x[2]^2/(2*key.scale^2))))/width^2
    }else{
      # define the error function in terms of pnorm
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # analytic expression for integral of g(x)/w when g(x) is half-normal
      int <- (1/width)*sqrt(pi/2)*key.scale*(-erf(x[1]/(key.scale*sqrt(2)))+
                                    erf(x[2]/(key.scale*sqrt(2))))
    }
    return(int)
  }else{
    return(integrate(distpdf, lower = x[1], upper = x[2], width=width,
                     ddfobj=ddfobj, select=select, index=index,
                     rel.tol = 1e-7,standardize=standardize,
                     stdint=stdint,point=point)$value)
  }
}
