#' Gradient of the integral of the detection fucntion (i.e., beta)
#' 
#' @param par.index the index of the parameter of interest
#' @param distance vector of distances
#' @param ddfobj the ddf object
#' @param width the truncation width
#' @param left the left truncation (defaults to zero)

detfct.grad <- function(par.index, distance, ddfobj, width = NULL, left = 0) {
  
  ## 1. Perform the preparations 
  ## ===========================
  # Key function
  key <- ddfobj$type
  
  pars <- getpar(ddfobj)
  par.indices <- getpar(ddfobj, index = TRUE)
  k <- sum(par.indices[-3])
  m <- par.indices[3] - k
  if (par.indices[2] != 0) key.scale <- pars[par.indices[2]] else key.scale <- NULL
  if (par.indices[1] != 0) key.shape <- pars[par.indices[1]] else key.shape <- NULL
  
  ## 2. Derive the gradients for parameters par.index > k
  ## ============================================
  zeros <- rep(0, length(distance))
  
  # If we are using adjustment terms.
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
  
  ## if the parameter is an adjustment parameter
  if (par.index > k) {
    par.index <- par.index - k
    # evaluate key function
    key.val <- switch(key,
                      hn    = keyfct.hn(distance, key.scale),
                      hr    = keyfct.hz(distance, key.scale, key.shape),
                      unif  = rep(1, length(distance)), # FTP: Why is this 1 and not 1/w?
                      gamma = keyfct.gamma(distance, key.scale, key.shape),
                      th1   = keyfct.th1(distance, key.scale, key.shape),
                      th2   = keyfct.th2(distance, key.scale, key.shape),
                      tpn   = keyfct.tpn(distance, ddfobj))
    
    
    ## Decide on adjustment term and run.
    adj.val <- switch(adj.series,
                      poly = adjfct.poly(distance, scaling, adj.order,
                                         adj.parm, adj.exp),
                      herm = adjfct.herm(distance, scaling, adj.order,
                                         adj.parm, adj.exp),
                      cos  = adjfct.cos(distance, scaling, adj.order,
                                        adj.parm, adj.exp))
    adj.term <- switch(adj.series,
                       poly = adj.poly(distance, scaling, adj.order[par.index] - 1),
                       herm = adj.herm(distance, scaling, adj.order[par.index] - 1),
                       cos  = adj.cos(distance, scaling, adj.order[par.index] - 1))
    
    ## Derive adj term and series, and key for distance = 0
    key.val.0 <- switch(key,
                        hn    = keyfct.hn(zeros, key.scale),
                        hr    = keyfct.hz(zeros, key.scale, key.shape),
                        unif  = 1,
                        gamma = keyfct.gamma(zeros, key.scale, key.shape),
                        th1   = keyfct.th1(zeros, key.scale, key.shape),
                        th2   = keyfct.th2(zeros, key.scale, key.shape))
    
    adj.val.0 <- switch(adj.series,
                        poly = adjfct.poly(zeros, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(zeros, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(zeros, scaling, adj.order,
                                          adj.parm, adj.exp))
    adj.term.0 <- switch(adj.series,
                         poly = adj.poly(zeros, scaling, adj.order[par.index] - 1),
                         herm = adj.herm(zeros, scaling, adj.order[par.index] - 1),
                         cos  = adj.cos(zeros, scaling, adj.order[par.index] - 1))
    
    # standardized value of the detection function
    grad <- key.val / (key.val.0 * (1 + adj.val.0)) ^ 2 * 
      (adj.term  * adj.val.0 - adj.val * adj.term.0)
  } else { # if a parameter is scale or shape
    ## when par is scale or shape, left side of the beta gradient is zero. we 
    ## only need to consider the right side.
    key.val.0 <- switch(key,
                        hn    = keyfct.hn(zeros, key.scale),
                        hr    = keyfct.hz(zeros, key.scale, key.shape),
                        unif  = 1, # FTP: Why is this 1 and not 1/w?
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
    
    D <- 1 / (key.val.0 * (1 + adj.val.0))
    
    if (par.index == 1) { # when parameter is scale
      scaled.dist.grad <- -(distance / key.scale ^ 2)
      
      key.val <- switch(key,
                        hn    = keyfct.hn(distance, key.scale),
                        hr    = keyfct.hz(distance, key.scale, key.shape),
                        unif  = 1)
      
      grad.series.adj.vals <- switch(
        adj.series,
        poly = grad.series.adj.poly(distance, key.scale, adj.order, adj.parm, adj.exp),
        herm = grad.series.adj.herm(distance, key.scale, adj.order, adj.parm, adj.exp),
        cos = grad.series.adj.cos(distance, key.scale, adj.order, adj.parm, adj.exp)
      )
      
      adj.val <- switch(adj.series,
                        poly = adjfct.poly(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(distance, scaling, adj.order,
                                          adj.parm, adj.exp))
      key.grad.val <- switch(key,
                             hn = keyfct.grad.hn(distance, key.scale),
                             hz = keyfct.grad.hz(distance, key.scale, key.shape))
      
      grad <- (key.val * grad.series.adj.vals * scaled.dist.grad + 
                 adj.val * key.grad.val) / (key.val.0 * (1 + adj.val.0))
      
    }
    else { # so when par.index = 2 and parameter is shape and key is hz
      # scaled.dist.grad <- 0
      # 
      key.grad.val <- keyfct.grad.hz(distance, key.scale, key.shape, shape = TRUE)
      
      adj.val <- switch(adj.series,
                        poly = adjfct.poly(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        herm = adjfct.herm(distance, scaling, adj.order,
                                           adj.parm, adj.exp),
                        cos  = adjfct.cos(distance, scaling, adj.order,
                                          adj.parm, adj.exp))
      
      grad <- adj.val * key.grad.val / (key.val.0 * (1 + adj.val.0))
    }
  }
  
  return(grad)
}
