#' Numerically integrate distance detection function over specified range
#' 
#' Computes integral of detection function over x for each observation.  The
#' method of computation depends on argument switches set and the type of
#' detection function.
#' 
#' If there are only intercepts in the scale model then the integrals can be
#' computed once and they apply to all observations.  If the scale model
#' includes variables other than intercepts then an integral must be computed
#' for each observation. Unless \code{doeachint} is set to TRUE, this is done
#' by creating a table of standardized integrals and interpolating the value
#' for each observation.  For the half-normal detection function the table is
#' computed at the very beginning and does not need to be recomputed.  For the
#' hazard rate the table needs to be recomputed for each iteration because the
#' shape function changes.  For small sample sizes it may be more efficient to
#' set doeachint=TRUE. For full independence models, integration of
#' p_1(x)*p_2(x) cannot use the table approach (except for the half-normal but
#' it has not been implemented) so each integral is computed. Except for this
#' latter integral, the integral can be computed by integrating a standardized
#' detection function and then multiplying by the scale for that observation.
#' This is a simple substitution approach for integration. y=x/scale.  However,
#' it does not work for the product because each has a separate scale; although
#' a common scale could be found algebraically at least for the half-normal
#' (see Laake (1999)).
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
integratedetfct <-
function(ddfobj,select,width,int.range,doeachint=FALSE,standardize=TRUE,point=FALSE)
{
#
# integratedetfct
#
# Computes integral of detection function over x for each observation.  The method of computation
# depends on switches set and the detection function.  If there are only inercepts in the scale model 
# then the integrals can be computed once and they apply to all observations.  If the scale model includes 
# variables other than intercepts then an integral must be computed for each observation. Unless doeachint
# is set to true, this is done by creating a table of standardized integrals and interpolating the value
# for each observation.  For the half-normal detection function the table is computed at the very beginning
# and does not need to be recomputed.  For the hazard rate the table needs to be recomputed for each iteration
# because the shape function changes.  For small sample sizes it may be more efficient to set doeachint=TRUE.
# For the FCI model, integration of g1(x)*g2(x) cannot use the table approach so each integral is computed.
#
# Arguments: 
# 
# ddfobj         - distance sampling object
# select         - selected values
# width          - scalar, vector or matrix of integration bounds
# int.range      - integration range
# doeachint      - logical that specifies whether each observation integral should be computed numerically
# standardize    - logical used to decide whether to divide through by the function evaluated at 0 
# point          - logical to determine if point count(TRUE) or line transect(FALSE)
#
# Value: 
#
# vector of integral values
#
#  Functions called:
#
#   integrate  - routine to integrate a function; used to integrate detfct or fcidupdetfc over x (0,width)          
#
  fpar=getpar(ddfobj)
#
#  Determine integration bounds
#
  if(is.vector(int.range)){
    if(length(int.range)==1){
      left <- 0
      right <- int.range
      samelimits <- TRUE
    }else if(length(int.range)==2){
      left <- int.range[1]
      right <- int.range[2]
      samelimits <- TRUE
    }else{
      stop("Invalid bounds setup")
    }
  
  }else{
    if(is.matrix(int.range)){
      samelimits <- FALSE
      left <- int.range[2:dim(int.range)[1],1]
      right <- int.range[2:dim(int.range)[1],2]
    }
  }
  if(ddfobj$intercept.only&samelimits)
  {
	  if(!point)
	  return(integrate(fx,lower=left,upper=right,ddfobj=ddfobj, select=c(TRUE,rep(FALSE,nrow(ddfobj$xmat)-1)),
                       width=width,standardize=standardize)$value)
    else
		return(integrate(fr,lower=left,upper=right,ddfobj=ddfobj, select=c(TRUE,rep(FALSE,nrow(ddfobj$xmat)-1)),
						width=width,standardize=standardize)$value)
  }
  else
    if(doeachint)
	{
		  integrals=rep(0,nrow(ddfobj$scale$dm[select,,drop=FALSE]))
		  xscale=scalevalue(ddfobj$scale$parameters, ddfobj$scale$dm[select,,drop=FALSE])
		  
		  for(i in 1:length(integrals))
		      if(!point)
			    integrals[i]=
					  xscale[i]*(gstdint((right/xscale)[i], xmin = 0, ddfobj=ddfobj, select=select, index=i, width=width,standardize=TRUE)-
					  gstdint((left/xscale)[i], xmin = 0, ddfobj=ddfobj, select=select, index=i, width=width,standardize=TRUE))
	          else
				  integrals[i]=
						  xscale[i]^2*(gstdint((right/xscale)[i], xmin = 0, ddfobj=ddfobj, select=select, index=i, width=width,standardize=TRUE,point=TRUE)-
						  gstdint((left/xscale)[i], xmin = 0, ddfobj=ddfobj, select=select, index=i, width=width,standardize=TRUE,point=TRUE))
		  return(integrals)
    }
	else
	{
      if(!is.null(ddfobj$shape))
        cgftab <- tablecgf(ddfobj,width=width,standardize=standardize,point)
      else
		cgftab=ddfobj$cgftab  
      # Calls predict.smooth.spline
		xscale=scalevalue(ddfobj$scale$parameters, ddfobj$scale$dm[select,,drop=FALSE])
        if(!point)
			  int=
                     xscale*(predict(cgftab, as.vector(right/xscale))$y -
			         predict(cgftab, as.vector(left/xscale))$y)
		else	
			  int=
			         xscale^2*(predict(cgftab, as.vector(right/xscale))$y  -
			         predict(cgftab, as.vector(left/xscale))$y) 
		return(int)
    }
}
