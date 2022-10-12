#'  Detection functions
#'
#' Various functions used to specify key and adjustment functions for
#' detection functions.
#'
#' Multi-covariate detection functions (MCDS) are represented by a function
#' \eqn{g(x,w,\theta)} where x is distance, z is a set of covariates and
#' \eqn{\theta} is the parameter vector.  The functions are defined such that
#' \eqn{g(0,w,\theta)=1} and the covariates modify the scale \eqn{(x/\sigma)}
#' where a log link is used to relate \eqn{\sigma} to the covariates,
#' \eqn{\sigma=exp(\theta*w)}. A CDS function is obtained with a constant
#' \eqn{\sigma} which is equivalent to an intercept  design matrix, z.
#'
#' \code{detfct} will call either a gamma, half-normal, hazard-rate or uniform
#' function only returning the probability of detection at that distance. In
#' addition to the simple model above, we may specify adjustment terms to fit
#' the data better. These adjustments are either Cosine, Hermite and simple
#' polynomials. These are specified as arguments to \code{detfct}, as detailed
#' below.
#'
#' \code{detfct} function which calls the others and assembles the final result
#' using either key(x)[1+series(x)] or
#' (key(x)[1+series(x)])/(key(0)[1+series(0)]) (depending on the value of
#' \code{standardize}).
#'
#' \code{keyfct.*} functions calculate key function values and \code{adjfct.*}
#' calculate adjustment term values.
#'
#' \code{scalevalue} for either detection function it computes the scale with
#' the log link using the parameters and the covariate design matrix
#'
#' \code{fx}, \code{fr} non-normalized probability density for line transects
#' and point counts respectively
#'
#' @aliases detfct adjfct.cos adjfct.herm hermite.poly adjfct.poly keyfct.hn
#'  keyfct.hz keyfct.gamma scalevalue fx fr distpdf
#' @usage detfct(distance, ddfobj, select=NULL, index=NULL, width=NULL,
#'               standardize = TRUE, stdint=FALSE, left=0)
#'
#' adjfct.cos(distance, scaling = 1, adj.order, adj.parm = NULL, adj.exp=FALSE)
#'
#' adjfct.poly(distance, scaling = 1, adj.order, adj.parm = NULL, adj.exp=FALSE)
#'
#' adjfct.herm(distance, scaling = 1, adj.order, adj.parm = NULL, adj.exp=FALSE)
#'
#' scalevalue(key.scale, z)
#'
#' keyfct.hn(distance, key.scale)
#'
#' keyfct.hz(distance, key.scale, key.shape)
#'
#' keyfct.gamma(distance, key.scale, key.shape)
#'
#' fx(distance,ddfobj,select=NULL,index=NULL,width=NULL,
#'    standardize=TRUE,stdint=FALSE, left=0)
#'
#' fr(distance,ddfobj,select=NULL,index=NULL,width=NULL,
#'    standardize=TRUE,stdint=FALSE)
#'
#' distpdf(distance,ddfobj,select=NULL,index=NULL,width=NULL,standardize=TRUE,
#'            stdint=FALSE,point=FALSE, left=0)
#'
#' @param distance  vector of distances
#' @param ddfobj distance sampling object (see \code{\link{create.ddfobj}})
#' @param z design matrix for scale function
#' @param select logical vector for selection of data values
#' @param index specific data row index
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#' @param adj.order vector of adjustment orders
#' @param adj.parm vector of adjustment parameters
#' @param width (right) truncation width
#' @param left (left) truncation distance
#' @param standardize logical used to decide whether to divide through by the
#' function evaluated at 0
#' @param scaling the scaling for the adjustment terms
#' @param stdint logical used to decide whether integral is standardized
#' @param point if TRUE, point counts; otherwise line transects
#' @param adj.exp if TRUE uses exp(adj) for adjustment to keep f(x)>0
#' @return
#' For \code{detfct}, the value is a vector of detection probabilities
#' For \code{keyfct.*}, vector of key function evaluations
#' For \code{adjfct.*}, vector of adjustment series evaluations
#' For \code{scalevalue}, vector of the scale parameters.
#' @author Jeff Laake, David L Miller
#' @seealso  \code{\link{mcds}},  \code{\link{cds}}
#' @references
#' Marques, F. F. C., & Buckland, S. T. (2003). Incorporating covariates into
#' standard line transect analyses. Biometrics, 59(4), 924-935.
#'
#' Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers, D.
#' L., & Thomas, L. (2004). Advanced Distance Sampling. Oxford University
#' Press, Oxford, UK.
#'
#' Becker, E. F. and P. X. Quang. 2009. A gamma-shaped detection function for
#' line transect surveys with mark-recapture and covariate data. Journal of
#' Agricultural Biological and Environmental Statistics 14:207-223.
#' @export detfct
#' @keywords internal
distpdf <- function(distance, ddfobj, select=NULL, index=NULL, width=NULL,
                    standardize=TRUE, stdint=FALSE, point=FALSE, left=0){

 if(!point){
   return(fx(distance=distance, ddfobj=ddfobj, select=select, index=index,
             width=width, standardize=standardize, stdint=stdint, left=left))
 }else
   return(fr(distance=distance, ddfobj=ddfobj, select=select, index=index,
             width=width, standardize=standardize, stdint=stdint))
}

fx <- function(distance, ddfobj, select=NULL, index=NULL, width=NULL,
               standardize=TRUE, stdint=FALSE, left=0){
  return(detfct(distance, ddfobj, select, index, width, standardize, stdint)/
          (width-left))
}

fr <- function(distance, ddfobj, select=NULL, index=NULL, width=NULL,
               standardize=TRUE, stdint=FALSE){
  return(detfct(distance, ddfobj, select, index, width, standardize, stdint)*
         2*distance/width^2)
}

detfct <- function(distance, ddfobj, select=NULL, index=NULL, width=NULL,
                   standardize=TRUE, stdint=FALSE, left=0){

  # Set of observations for computation of detection function can
  # be specified with logical (select) and numeric (index) values.
  # Either or both can be specified although the latter is unlikely.
  if(is.null(select)){
    # use all
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm
      shape.dm <- ddfobj$shape$dm
    }else{
      # use only those with specific indices
      scale.dm <- ddfobj$scale$dm[index, , drop=FALSE]
      shape.dm <- ddfobj$shape$dm[index, , drop=FALSE]
    }
    ddfobj$scale$dm <- scale.dm
    ddfobj$shape$dm <- shape.dm
  }else{
    # Use those with select=TRUE
    if(is.null(index)){
      scale.dm <- ddfobj$scale$dm[select, , drop=FALSE]
      shape.dm <- ddfobj$shape$dm[select, , drop=FALSE]
    }else{
      # use the numeric index within those with select=TRUE
      scale.dm <- ddfobj$scale$dm[select, , drop=FALSE][index, , drop=FALSE]
      shape.dm <- ddfobj$shape$dm[select, , drop=FALSE][index, , drop=FALSE]
    }
    ddfobj$scale$dm <- scale.dm
    ddfobj$shape$dm <- shape.dm
  }

  # Key function
  key <- ddfobj$type

  # calculate the key scale
  if(stdint){
    if(is.null(index)){
      key.scale <- rep(1, nrow(scale.dm))
    }else{
      key.scale <- 1
    }
  }else{
    if(!is.null(ddfobj$scale)){
      key.scale <- scalevalue(ddfobj$scale$parameters, scale.dm)
    }
  }

  # calculate the key shape
  if(!is.null(ddfobj$shape)){
    key.shape <- scalevalue(ddfobj$shape$parameters, shape.dm)
  }

  # for gamma shape parameter must be >1, see Becker and Quang (2009) p 213
  if(key=="gamma"){
    key.shape <- key.shape + 1
    key.shape[key.shape==1] <- key.shape[key.shape==1] + 0.000001
  }

  # evaluate key function
  g <- switch(key,
              hn    = keyfct.hn(distance, key.scale),
              hr    = keyfct.hz(distance, key.scale, key.shape),
              unif  = rep(1, length(distance)),
              gamma = keyfct.gamma(distance, key.scale, key.shape),
              th1   = keyfct.th1(distance, key.scale, key.shape),
              th2   = keyfct.th2(distance, key.scale, key.shape),
              tpn   = keyfct.tpn(distance, ddfobj))

  # If we are using adjustment terms.
  if(!is.null(ddfobj$adjustment)){
    adj.series <- ddfobj$adjustment$series
    adj.scale <- ddfobj$adjustment$scale
    adj.order <- ddfobj$adjustment$order
    adj.parm <- ddfobj$adjustment$parameters
    adj.exp <- ddfobj$adjustment$exp

    # Find out if we are scaling by width or by key scale
    if(adj.scale == "width"){
      scaling <- width
    }else{
      scaling <- key.scale
    }

    ## Decide on adjustment term and run.
    adj.vals <- switch(adj.series,
                       poly = adjfct.poly(distance, scaling, adj.order,
                                          adj.parm, adj.exp),
                       herm = adjfct.herm(distance, scaling, adj.order,
                                          adj.parm, adj.exp),
                       cos  = adjfct.cos(distance, scaling, adj.order,
                                         adj.parm, adj.exp))

    # calculate detection function with adjustments
    g <- g * (1 + adj.vals)

    # If we have adjustment terms then we need to standardize the detection
    # function. So find the values for the key and adjustment terms at 0
    # this cancels in the likelihood, so we don't need it in optimisation
    if(standardize){
      if(key == "gamma"){
        # for the gamma, use apex.gamma to find the apex first, then eval
        # need to update the scale to be +1 in this apex call
        ddfobj$shape$parameters <- log(exp(ddfobj$shape$parameters)+1)
        zeros <- as.vector(apex.gamma(ddfobj))[1]
      }else{
        zeros <- rep(0, length(distance))
      }

      key.val.0 <- switch(key,
                          hn    = keyfct.hn(zeros, key.scale),
                          hr    = keyfct.hz(zeros, key.scale, key.shape),
                          unif  = rep(1, length(zeros)),
                          gamma = keyfct.gamma(zeros, key.scale, key.shape),
                          th1   = keyfct.th1(zeros, key.scale, key.shape),
                          th2   = keyfct.th2(zeros, key.scale, key.shape))

      # now compute adjustments
      zeros <- rep(0, length(distance))

      adj.val.0 <- switch(adj.series,
                          poly = adjfct.poly(zeros, scaling, adj.order,
                                             adj.parm, adj.exp),
                          herm = adjfct.herm(zeros, scaling, adj.order,
                                             adj.parm, adj.exp),
                          cos  = adjfct.cos(zeros, scaling, adj.order,
                                            adj.parm, adj.exp))

      # standardized value of the detection function
      g <- g/(key.val.0 * (1 + adj.val.0))
    } # end standardize
  } # end adjustments
  return(g)
}
