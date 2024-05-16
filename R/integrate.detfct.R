#' Integrates the (non-)standardised detection function. The output is beta in 
#' the documentation.
#' 
#' look at int1 <- integratepdf(ddfobj, select=!x$binned, width=right,
#' int.range=intrange[!x$binned, ],
#' doeachint=doeachint,
#' point=misc.options$point, standardize=standardize,
#' left=left)
#' in flpt.lnl.r for inspiration on how to use this function. 

integrate.detfct <- function(int.range, ddfobj, width, standardize, 
                 point = FALSE, left = 0, stdint = FALSE, select = NULL,
                 index = 1) {
  ## The integral over distance of the non-standardised detection function
  g <- integrate(detfct, lower = int.range[1], upper = int.range[2],
                 ddfobj = ddfobj, width = width, standardize = standardize,
                 select = select, stdint = stdint, left = left)$value
  
  ## Version below is copied from previous code in mrds. 
  # B <- gstdint(int.range[1, ], ddfobj = ddfobj, index = index, select = select,
  #              width = width, standardize = standardize, point = point,
  #              stdint = stdint, left = left)
  # 
  # detfct(1, ddfobj = ddfobj, width = width, standardize = standardize,
  #        select = select, stdint = stdint, left = left, index = 1)
  
  ## If the detection function is standardised
  if (standardize) {
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
    

    key.val.0 <- switch(key,
                        hn    = keyfct.hn(0, key.scale),
                        hr    = keyfct.hz(0, key.scale, key.shape),
                        unif  = 1,
                        gamma = keyfct.gamma(0, key.scale, key.shape),
                        th1   = keyfct.th1(0, key.scale, key.shape),
                        th2   = keyfct.th2(0, key.scale, key.shape))
    
    adj.val.0 <- switch(adj.series,
                        poly = adjfct.poly(0, scaling, adj.order, adj.parm, 
                                           adj.exp),
                        herm = adjfct.herm(0, scaling, adj.order, adj.parm, 
                                           adj.exp),
                        cos  = adjfct.cos(0, scaling, adj.order, adj.parm, 
                                          adj.exp))
    
    return(g / (key.val.0 * (1 + adj.val.0)))
    
  } else {
    return(g)
  }
}
