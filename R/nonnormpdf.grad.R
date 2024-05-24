#' Gradient of the detection function
#'
#' Various functions used to specify key and adjustment functions for
#' gradients of detection functions.
#' 
#' So far, only developed for the half-normal, hazard-rate and uniform key
#' functions in combination with cosine, simple polynomial and Hermite 
#' polynomial adjustments. It is only called by the gradient-based solver
#' and not accessible to the general user
#'
#' \code{detfct.grad} will call either a half-normal, hazard-rate or uniform
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
#' Defaults to TRUE.
#' 
#'@section Dependencies: 
#'  \describe {
#'    \item{keyfct.XX}
#'    \item{keyfct.grad.XX}
#'    \item{adj.YY}
#'    \item{adjfct.YY}
#'    \item{adj.grd.YY}
#'  }
nonnormpdf.grad <- function(distance, par.index, ddfobj, standardize = FALSE, 
                            point = FALSE, width = NULL, left = 0) {
  
  ## 1. Extract information from ddfobj and prepare for gradient evaluation
  ## ======================================================================
  # Key function
  key <- ddfobj$type
  
  pars <- getpar(ddfobj)
  par.indices <- getpar(ddfobj, index = TRUE)
  k <- sum(par.indices[-3])
  m <- par.indices[3] - k
  
  if(par.indices[2] != 0) {
    key.scale <- pars[par.indices[2]]
  } else {
    key.scale <- NULL
  }
  if(par.indices[1] != 0) {
    key.shape <- pars[par.indices[1]]
  } else {
    key.shape <- NULL
  }
  # key.scale <- ifelse(par.indices[2] != 0, pars[par.indices[2]], NULL)
  # key.shape <- ifelse(par.indices[1] != 0, pars[par.indices[1]], NULL)
  
  
  ## 2. Derive the gradients for parameter if index > k, i.e., it is an 
  ##    adjustment parameter.
  ## ==================================================================
  # zeros <- rep(0, length(distance)) # is distance not always length 1? 
                                      # I think so, so commented out
  
  ## Extact the information about the adjustment term
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
  

  ## Evaluate the gradient of the standardised detection function if 
  ## standardize == FALSE. Things are a lot easier then.
  if (!standardize) { # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< standardize == FALSE
    ## IF LINE TRANSECT DATA
    if (!point) {
      ## If the parameter is an adjustment parameter
      if (par.index > k) {
        ## Derive the adjustment parameter index (j' in the documentation)
        adj.par.index <- par.index - k
        
        ## Evaluate the specified key function
        ## Currently implemented for: half-normal, hazard-rate, uniform
        key.val <- switch(key,
                          hn    = keyfct.hn(distance, key.scale),
                          hr    = keyfct.hz(distance, key.scale, key.shape),
                          unif  = 1, 
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
        grad <- key.val * adj.term 
        
      } 
      else {
        ## If the parameter is a key function parameter (i.e., scale or shape)
        ## Evaluate the key function and adjustment series for distance = 0
        if (par.index == 1) { ## If the parameter is the scale parameter
          ## Derive the gradient of the scaled distance w.r.t. the parameter
          scaled.dist.grad <- -(distance / key.scale ^ 2)
          
          ## Evaluate the key function and adjustment series
          key.val <- switch(key,
                            hn    = keyfct.hn(distance, key.scale),
                            hr    = keyfct.hz(distance, key.scale, key.shape),
                            unif  = 1, 
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
                                 hz = keyfct.grad.hz(distance, key.scale, 
                                                     key.shape),
                                 unif = 0)
          
          ## Derive the gradient of the non-standardised detection function 
          grad <- (key.val * grad.adj.series.val * scaled.dist.grad + 
            adj.val * key.grad.val)
          
        }
        else { ## If par.index = 2, the parameter is shape and key is hazard-rate.
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
          grad <- adj.val * key.grad.val 
        }
      }
    } 
    ## IF POINT TRANSECT DATA
    else {
      ## If the parameter is an adjustment parameter
      if (par.index > k) {
        ## Derive the adjustment parameter index (j' in the documentation)
        adj.par.index <- par.index - k
        
        ## Evaluate the specified key function
        ## Currently implemented for: half-normal, hazard-rate, uniform
        key.val <- switch(key,
                          hn    = keyfct.hn(distance, key.scale),
                          hr    = keyfct.hz(distance, key.scale, key.shape),
                          unif  = 1, 
                          gamma = keyfct.gamma(distance, key.scale, key.shape),
                          th1   = keyfct.th1(distance, key.scale, key.shape),
                          th2   = keyfct.th2(distance, key.scale, key.shape),
                          tpn   = keyfct.tpn(distance, ddfobj))
        
        
        ## Evaluate the specified adjustment term 
        adj.term <- switch(adj.series,
                           poly = adj.poly(distance, scaling, 
                                           adj.order[adj.par.index] - 1),
                           herm = adj.herm(distance, scaling, 
                                           adj.order[adj.par.index] - 1),
                           cos  = adj.cos(distance, scaling, 
                                          adj.order[adj.par.index] - 1))
        
        ## Derive the gradient of the non-standardised detection function
        grad <- key.val * adj.term 
        
      } else { 
        ## If the parameter is a key function parameter (i.e., scale or shape)
        ## Evaluate the key function and adjustment series for distance = 0
        
        if (par.index == 1) { ## If the parameter is the scale parameter
          ## Derive the gradient of the scaled distance w.r.t. the parameter
          scaled.dist.grad <- -(distance / key.scale ^ 2)
          
          ## Evaluate the key function and adjustment series
          key.val <- switch(key,
                            hn    = keyfct.hn(distance, key.scale),
                            hr    = keyfct.hz(distance, key.scale, key.shape),
                            unif  = 1, 
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
                                 hz = keyfct.grad.hz(distance, key.scale, 
                                                     key.shape),
                                 unif = 0)
          
          ## Derive the gradient of the non-standardised detection function 
          grad <- distance * key.val * grad.adj.series.val * scaled.dist.grad + 
            adj.val * key.grad.val
          
        }
        else { ## If par.index = 2, the parameter is shape and key is hazard-rate.
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
          grad <- distance * adj.val * key.grad.val 
        }
      }
    }
    
  }   
  ## Evaluate the gradient of the standardised detection function if 
  ## standardize == TRUE
  else { # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< standardize == TRUE
    ## If the parameter is an adjustment parameter
    if (par.index > k) {
      ## Derive the adjustment parameter index (j' in the documentation)
      adj.par.index <- par.index - k
      
      ## Evaluate the specified key function
      ## Currently implemented for: half-normal, hazard-rate, uniform
      key.val <- switch(key,
                        hn    = keyfct.hn(distance, key.scale),
                        hr    = keyfct.hz(distance, key.scale, key.shape),
                        unif  = 1, 
                        gamma = keyfct.gamma(distance, key.scale, key.shape),
                        th1   = keyfct.th1(distance, key.scale, key.shape),
                        th2   = keyfct.th2(distance, key.scale, key.shape),
                        tpn   = keyfct.tpn(distance, ddfobj))
      
      
      ## Evaluate the specified adjustment series and term 
      adj.val <- switch(adj.series,
                        poly = adjfct.poly(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(distance, scaling, adj.order,
                                          adj.parm, adj.exp))
      adj.term <- switch(adj.series,
                         poly = adj.poly(distance, scaling, 
                                         adj.order[adj.par.index] - 1),
                         herm = adj.herm(distance, scaling, 
                                         adj.order[adj.par.index] - 1),
                         cos  = adj.cos(distance, scaling, 
                                        adj.order[adj.par.index] - 1))
      
      ## Evaluate the specified key function for distance = 0
      key.val.0 <- switch(key,
                          hn    = keyfct.hn(0, key.scale),
                          hr    = keyfct.hz(0, key.scale, key.shape),
                          unif  = 1, 
                          gamma = keyfct.gamma(0, key.scale, key.shape),
                          th1   = keyfct.th1(0, key.scale, key.shape),
                          th2   = keyfct.th2(0, key.scale, key.shape))
      
      ## Evaluate the specified adjustment series and term for distance = 0
      adj.val.0 <- switch(adj.series,
                          poly = adjfct.poly(0, scaling, adj.order,
                                             adj.parm, adj.exp),
                          herm = adjfct.herm(0, scaling, adj.order,
                                             adj.parm, adj.exp),
                          cos  = adjfct.cos(0, scaling, adj.order,
                                            adj.parm, adj.exp))
      adj.term.0 <- switch(adj.series,
                           poly = adj.poly(0, scaling, 
                                           adj.order[adj.par.index] - 1),
                           herm = adj.herm(0, scaling, 
                                           adj.order[adj.par.index] - 1),
                           cos  = adj.cos(0, scaling, 
                                          adj.order[adj.par.index] - 1))
      
      ## Derive the gradient of the standardised detection function 
      grad <- key.val / (key.val.0 * (1 + adj.val.0)) ^ 2 * 
        (adj.term  * (1 + adj.val.0) - adj.val * (1 + adj.term.0)) ## CHECK THIS!
    } else { 
      ## If the parameter is a key function parameter (i.e., scale or shape)
      ## Evaluate the key function and adjustment series for distance = 0
      key.val.0 <- switch(key,
                          hn    = keyfct.hn(0, key.scale),
                          hr    = keyfct.hz(0, key.scale, key.shape),
                          unif  = 1, 
                          gamma = keyfct.gamma(0, key.scale, key.shape),
                          th1   = keyfct.th1(0, key.scale, key.shape),
                          th2   = keyfct.th2(0, key.scale, key.shape),
                          tpn   = keyfct.tpn(0, ddfobj))
      adj.val.0 <- switch(adj.series,
                          poly = adjfct.poly(0, scaling, adj.order,
                                             adj.parm, adj.exp),
                          herm = adjfct.herm(0, scaling, adj.order,
                                             adj.parm, adj.exp),
                          cos  = adjfct.cos(0, scaling, adj.order,
                                            adj.parm, adj.exp))
      
      if (par.index == 1) { ## If the parameter is the scale parameter
        ## Derive the gradient of the scaled distance w.r.t. the parameter
        scaled.dist.grad <- -(distance / key.scale ^ 2)
        
        ## Evaluate the key function and adjustment series
        key.val <- switch(key,
                          hn    = keyfct.hn(distance, key.scale),
                          hr    = keyfct.hz(distance, key.scale, key.shape),
                          unif  = 1, 
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
        
        ## Evaluate the gradient of the key function w.r.t. the parameter
        key.grad.val <- switch(key,
                               hn = keyfct.grad.hn(distance, key.scale),
                               hz = keyfct.grad.hz(distance, key.scale, 
                                                   key.shape),
                               unif = 0)
        
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
        
        ## Derive the gradient of the standardised detection function 
        grad <- (key.val * grad.adj.series.val * scaled.dist.grad + 
                   adj.val * key.grad.val) / (key.val.0 * (1 + adj.val.0))
        
      }
      else { ## If par.index = 2, the parameter is shape and key is hazard-rate.
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
        
        ## Derive the gradient of the standardised detection function 
        grad <- adj.val * key.grad.val / (key.val.0 * (1 + adj.val.0))
      }
    }
  }
  
  ## Return the gradient
  return(grad)
}
