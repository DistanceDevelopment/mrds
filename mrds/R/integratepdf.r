#' Numerically integrate pdf of observed distances over specified ranges
#' 
#' Computes integral of pdf of observed distances over x for each observation.  The
#' method of computation depends on argument switches set and the type of
#' detection function.
#' 
#' If either doeachint is set or there is only one integral then they are computed 
#' using integrate; otherwise, it uses the cgftab which is a spline
#' fitted to a table of standardized integrals and the value is interpolated from the
#' spline for each observation.
#'
#' @param ddfobj distance detection function specification
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param int.range integration range
#' @param doeachint logical that specifies whether each observation integral
#'   should be computed numerically
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count(TRUE) or line
#'   transect(FALSE)
#' @return vector of integral values - one for each observation
#' @author Jeff Laake
#' @keywords utility
integratepdf <- function(ddfobj, select, width, int.range, doeachint=FALSE,
                         standardize=TRUE, point=FALSE){
#  Make sure there is consistency between integration ranges and data
#  It is ok to have a single observation with multiple ranges or a single range
#  with multiple observations but otherwise the numbers must agree if both >1

  if(is.vector(int.range)){
    stop("\nInternal error - int.range not a matrix")
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
  # Now either compute integral by integrating function or by use of spline in
  # cgftab.  If either doeachint=TRUE or only a single integral, then integrate
  # function using gstdint with stdint=FALSE
  if(doeachint || (nobs==1 & nrow(int.range)==1)){
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
                               select=NULL,width=width, standardize=standardize,
                               point=point, stdint=FALSE)
          }
        }else{
          integrals[i] <- gstdint(int.range[i,], ddfobj=ddfobj, index=index[1],
                                  select=NULL, width=width,
                             standardize=standardize, point=point, stdint=FALSE)
        }
      }
    }
    return(integrals)
  }else{
    # otherwise, use cgftab spline.
    cgftab <- ddfobj$cgftab
    if(!is.null(ddfobj$scale)){
      xscale <- scalevalue(ddfobj$scale$parameters,
                           ddfobj$scale$dm[select,,drop=FALSE])
    }else{
      xscale <- 1
    }

    if(!point){
      integrals <- xscale*(predict(cgftab, as.vector(int.range[,2]/xscale))$y -
                           predict(cgftab, as.vector(int.range[,1]/xscale))$y)
    }else{
      integrals <- xscale^2*(predict(cgftab,as.vector(int.range[,2]/xscale))$y -
                             predict(cgftab, as.vector(int.range[,1]/xscale))$y)
    }
    return(integrals)
  }
}
