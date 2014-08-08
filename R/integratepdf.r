#' Numerically integrate pdf of observed distances over specified ranges
#'
#' Computes integral of pdf of observed distances over x for each observation.
#' The method of computation depends on argument switches set and the type of
#' detection function.
#'
#' @param ddfobj distance detection function specification
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param int.range integration range matrix; vector is converted to matrix
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count (\code{TRUE}) or line
#'   transect (\code{FALSE})
#' @return vector of integral values - one for each observation
#' @author Jeff Laake
#' @keywords utility
integratepdf <- function(ddfobj, select, width, int.range,
                         standardize=TRUE, point=FALSE){
  # Make sure there is consistency between integration ranges and data
  # It is ok to have a single observation with multiple ranges or a single range
  # with multiple observations but otherwise the numbers must agree if both >1

  if(!is.matrix(int.range)){
    if(is.vector(int.range) && length(int.range)==2){
      int.range=matrix(int.range,ncol=2,nrow=1)
    }else{
        stop("\nInternal error - int.range is not a matrix and cannot be converted to the required matrix structure")
    }
  }

  if(is.null(select)){
    nobs <- nrow(ddfobj$xmat)
    index <- 1:nobs
  }else{
    nobs <- sum(select)
    index <- which(select)
  }

  if(nrow(int.range)>1 && nobs>1 && nrow(int.range)!=nobs){
    stop("\n Number of integration ranges (int.range) does not match number of observations\n")
  }

  # Now compute integrals
  if(nobs==1){
    return(gstdint(int.range[1,], ddfobj=ddfobj, index=1, select=NULL,
           width=width, standardize=standardize, point=point, stdint=FALSE))
  }else{
    nintegrals <- max(nobs,nrow(int.range))
    integrals <- vector("numeric",nintegrals)
    for(i in 1:nintegrals){
      if(nobs>1){
        if(nrow(int.range)>1){
          integrals[i] <- gstdint(int.range[i,], ddfobj=ddfobj,
                                  index=index[i], select=NULL, width=width,
                                  standardize=standardize, point=point,
                                  stdint=FALSE)
        }else{
          integrals[i] <- gstdint(int.range, ddfobj=ddfobj, index=index[i],
                                  select=NULL,width=width,
                                  standardize=standardize,point=point,
                                  stdint=FALSE)
        }
      }else{
        integrals[i] <- gstdint(int.range[i,], ddfobj=ddfobj, index=index[1],
                                select=NULL, width=width,
                                standardize=standardize, point=point,
                                stdint=FALSE)
      }
    }
  }
  return(integrals)
}
