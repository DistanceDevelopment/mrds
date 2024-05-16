#' Gradient of the integral of the detection function, i.e., d beta/d theta in 
#' the documentation. This gradient of the integral is the same as the integral
#' of the gradient, thanks to Leibniz integral rule. 
#' 
#' @param par.index the index of the parameter of interest
#' @param distance vector of distances
#' @param ddfobj the ddf object
#' @param width the truncation width
#' @param left the left truncation (defaults to zero)

integrate.detfct.grad <- function(distance, par.index, ddfobj, int.range, width, 
                      standardize, point = FALSE, left = 0, 
                      stdint = FALSE, select = NULL, index = NULL) {
  
  ## If the non-standardised detection function is required, simply 
  ## integrate the gradient of the detection function directly.
  if (!standardize) {
    out <- integrate(detfct.grad, lower = int.range[1], upper = int.range[2],
                     par.index = par.index, ddfobj = ddfobj, width = width,
                     standardize = FALSE)$value
  } 
  ## If the standardised detection function is specified, things get a little
  ## more complicated.
  else {
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
    
    ## If we are using adjustment terms.
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
    
    ## Create zero vector of same length as distance
    zeros <- rep(0, length(distance))
    
    ## Evaluate the key function and adjustment series for distance = 0
    key.val.0 <- switch(key,
                        hn    = keyfct.hn(zeros, key.scale),
                        hr    = keyfct.hz(zeros, key.scale, key.shape),
                        unif  = rep(1, length(distance)), 
                        gamma = keyfct.gamma(zeros, key.scale, key.shape),
                        th1   = keyfct.th1(zeros, key.scale, key.shape),
                        th2   = keyfct.th2(zeros, key.scale, key.shape),
                        tpn   = keyfct.tpn(zeros, ddfobj))
    adj.val.0 <- switch(adj.series,
                        poly = adjfct.poly(zeros, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(zeros, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(zeros, scaling, adj.order,
                                          adj.parm, adj.exp))
    
    ## We evaluate the four parts. These are referred to as part A through D. 
    ## Part D: 1 over the detection function at distance 0. 
    ## ======================================================================
    D <- 1 / (key.val.0 * (1 + adj.val.0))
    
    ## Part C: integrate the gradient of the detection function
    ## ========================================================
    C <- integrate(detfct.grad, lower = int.range[1], upper = int.range[2],
                   par.index = par.index, ddfobj = ddfobj, width = width,
                   standardize = FALSE)$value 
    # Note: Cis the same regardless of standardization.. is that correct?
    
    ## 2. Derive the gradients for parameters par.index > k
    ## ====================================================
    ## if the parameter is an adjustment parameter
    if (par.index > k) {
      
      ## Part A: the derivative of the series w.r.t. an adjustment parameter
      ## ===================================================================
      A1 <- - 1 / (key.val.0 * (adj.val.0) ^ 2)
      A2 <- switch(adj.series, 
                   poly = 0,
                   herm = hermite.poly(0, adj.order[par.index - k]),
                   cos = 1)
      A <- A1 * A2
      
      ## Part B: the integral over distance of the non-standardised detection 
      ##         function
      ## =====================================================================
      B <- integrate(detfct, lower = int.range[1], upper = int.range[2],
                     ddfobj = ddfobj, width = width, standardize = standardize,
                     select = select, stdint = stdint, left = left)$value
      
      ## Version below is copied from previous code in mrds. 
      # B <- gstdint(int.range[1, ], ddfobj = ddfobj, index = index, select = select,
      #              width = width, standardize = standardize, point = point,
      #              stdint = stdint, left = left)

    } else { # if a parameter is scale or shape
      
      ## Part A*B: the derivative of the series w.r.t a key parameter is zero 
      ##           when distance is zero. Therefore part A*B becomes zero.
      ## ========================================================================
      A <- B <- 0
    }
    out <- A * B + C * D
  }
  
  return(out)
}
