#' Gradient of the non-normalised pdf of distances or the detection function
#' for the distances.
#' 
#' This function has been updated to match distpdf closely, so
#' that it has the same flexibility. Effectively, it gives the gradient of 
#' distpdf or detfct, whichever one is specified. 
#' 
#' Various functions used to specify key and adjustment functions for
#' gradients of detection functions.
#' 
#' So far, only developed for the half-normal, hazard-rate and uniform key
#' functions in combination with cosine, simple polynomial and Hermite 
#' polynomial adjustments. It is only called by the gradient-based solver
#' and should not be called by the general user.
#'
#' \code{distpdf.grad} will call either a half-normal, hazard-rate or uniform
#' function with adjustment terms to fit the data better, returning the 
#' gradient of detection at that distance w.r.t. the parameters. The adjustments 
#' are either cosine, Hermite or simple polynomial.
#' 
#' @param par.index the index of the parameter of interest
#' @param distance vector of distances
#' @param ddfobj the ddf object
#' @param width the truncation width
#' @param left the left truncation (defaults to zero)
#' @param standardize whether the function should return the gradient of the  
#' standardized detection function g(x)/g(0) (TRUE), or simply of g(0) (FALSE). 
#' Currently only implemented for standardize = FALSE. 
#' 
#' @return the gradient of the non-normalised pdf or detection w.r.t. to 
#' the parameter with parameter index \code{par.index}. 
#' 
#'@section Dependencies: 
#'  \describe {
#'    \item{keyfct.XX}
#'    \item{keyfct.grad.XX}
#'    \item{adj.YY}
#'    \item{adjfct.YY}
#'    \item{adj.grd.YY}
#'  }
#' 
#' @author Felix Petersma
distpdf.grad <- function(distance, par.index, ddfobj, standardize = FALSE, 
                        width, point, left = 0, pdf.based = TRUE) {
  
  ## 1. Extract information from ddfobj and prepare for gradient evaluation
  ## ======================================================================
  ## Key function
  key <- ddfobj$type

  pars <- getpar(ddfobj)
  par.indices <- getpar(ddfobj, index = TRUE)
  k <- sum(par.indices[-3] > 0)
  m <- par.indices[3] - k
  
  if(par.indices[2] != 0) {
    key.scale <- scalevalue(pars[par.indices[2]], ddfobj$scale$dm[1])
  } else {
    key.scale <- NULL
  }
  if(par.indices[1] != 0) {
    key.shape <- scalevalue(pars[par.indices[1]], ddfobj$shape$dm[1])
  } else {
    key.shape <- NULL
  }
  
  ## determine which parameter it is
  if (par.index > k) {
    par.type <- "adjustment"
  } else if (!is.null(key.shape) & par.index == 1) {
    par.type <- "shape"
  } else {
    par.type <- "scale"
  }

  zeros <- 0 # rep(0, length(distance)) 
  
  ## Extract the information about the adjustment term
  adj.series <- ddfobj$adjustment$series
  adj.scale <- ddfobj$adjustment$scale
  adj.order <- ddfobj$adjustment$order
  adj.parm <- ddfobj$adjustment$parameters
  adj.exp <- ddfobj$adjustment$exp
  
  ## Find out if we are scaling by width or by key scale
  if(adj.scale == "width"){
    scaling <- width
  }else{
    scaling <- key.scale
  }
  
  ## If the parameter is an adjustment parameter
  if (par.type == "adjustment") {
    ## Derive the adjustment parameter index (j' in the documentation)
    adj.par.index <- par.index - k
    
    ## Evaluate the specified key function
    ## Currently implemented for: half-normal, hazard-rate, uniform
    key.val <- switch(key,
                      hn    = keyfct.hn(distance, key.scale),
                      hr    = keyfct.hz(distance, key.scale, key.shape),
                      unif  = 1, # rep(1, length(distance)), 
                      gamma = keyfct.gamma(distance, key.scale, key.shape),
                      th1   = keyfct.th1(distance, key.scale, key.shape),
                      th2   = keyfct.th2(distance, key.scale, key.shape),
                      tpn   = keyfct.tpn(distance, ddfobj))
    
    
    ## Evaluate the specified adjustment term 
    adj.term <- switch(adj.series,
                       poly = adj.poly(distance, scaling, 
                                       adj.order[adj.par.index]),
                       herm = adj.herm(distance, scaling, 
                                       adj.order[adj.par.index]),
                       cos  = adj.cos(distance, scaling, 
                                      adj.order[adj.par.index]))
    
    ## Derive the gradient of the non-standardised detection function
    if (point) {
      grad <- key.val * adj.term
      if (pdf.based) grad <- grad * 2 * distance / (width ^ 2)
    } else {
      grad <- key.val * adj.term
      if (pdf.based) grad <- grad  / (width - left)
    }
    
  } else { 
    ## If the parameter is a key function parameter (i.e., scale or shape)
    ## Evaluate the key function and adjustment series for distance = 0
    
    if (par.type == "scale") { ## If the parameter is the scale parameter
      ## Derive the gradient of the scaled distance w.r.t. the parameter
      ## This is zero when scaling is done using width instead of the scale par
      if (adj.scale == "scale") {
        scaled.dist.grad <- -1 * distance / (key.scale ^ 2)
      } else{
        scaled.dist.grad <- 0
      }
            
      ## Evaluate the key function and adjustment series
      key.val <- switch(key,
                        hn    = keyfct.hn(distance, key.scale),
                        hr    = keyfct.hz(distance, key.scale, key.shape),
                        unif  = 1, #rep(1, length(distance)), 
                        gamma = keyfct.gamma(distance, key.scale, key.shape),
                        th1   = keyfct.th1(distance, key.scale, key.shape),
                        th2   = keyfct.th2(distance, key.scale, key.shape),
                        tpn   = keyfct.tpn(distance, ddfobj))
      adj.val <- switch(adj.series,
                        poly = adjfct.poly(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(distance, scaling, adj.order,
                                          adj.parm, adj.exp))
      
      ## Evaluate gradient of the adjustment series w.r.t. the scaled distance
      grad.adj.series.val <- switch(
        adj.series,
        poly = grad.adj.series.poly(distance, key.scale, adj.order, adj.parm, 
                                    adj.exp),
        herm = grad.adj.series.herm(distance, key.scale, adj.order, adj.parm, 
                                    adj.exp),
        cos = grad.adj.series.cos(distance, key.scale, adj.order, adj.parm, 
                                  adj.exp)
      )
      
      ## Evaluate the gradient of the key function w.r.t. the parameter
      key.grad.val <- switch(key,
                             hn = keyfct.grad.hn(distance, key.scale),
                             hr = keyfct.grad.hz(distance, key.scale, 
                                                 key.shape),
                             unif = 0) # rep(0, distance))
      
      ## Derive the gradient of the non-standardised detection function
      if (point) {
        grad <- (key.val * grad.adj.series.val * scaled.dist.grad + 
                   key.grad.val * (1 + adj.val)) * key.scale
        if (pdf.based) grad <- grad * 2 * distance / (width ^ 2)
      } else {
        grad <- (key.val * grad.adj.series.val * scaled.dist.grad + 
                   key.grad.val * (1 + adj.val)) * key.scale
        if (pdf.based) grad <- grad / (width - left)
      }
    }
    else { ## If par.type == "shape", the parameter is shape and key is hazard-rate.
      ## Since the derivative of the scaled distance wrt the shape 
      ## parameter is 0, the majority of the equation becomes zero.
      
      ## Evaluate the derivative of the key function w.r.t. the shape parameter
      key.grad.val <- keyfct.grad.hz(distance, key.scale, key.shape, 
                                     shape = TRUE)
      
      ## Evaluate the adjustment series
      adj.val <- switch(adj.series,
                        poly = adjfct.poly(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(distance, scaling, adj.order,
                                          adj.parm, adj.exp))
      
      ## Derive the gradient of the non-standardised detection function 
      if (point) {
        grad <- (key.grad.val * (1 + adj.val)) * key.shape
        if (pdf.based) grad <- grad * 2 * distance / (width ^ 2)
        
      } else {
        grad <- (key.grad.val * (1 + adj.val)) * key.shape
        if (pdf.based) grad <- grad / (width - left)
      }
    }
  }
  
  ## Return the gradient
  return(grad)
}

# ## A wrapper in case it is important that the gradient is positive
# pdistpdf.grad <- function(distance, par.index, ddfobj, standardize, 
#                           width, left = 0) {
#   res <- distpdf.grad(distance = distance, par.index = par.index, ddfobj=ddfobj, 
#                       standardize=standardize, width = width, #point=point, # point not yet implemented
#                       left = left)
#   res[res < 1e-16] <- 0
#   return(res)
# }