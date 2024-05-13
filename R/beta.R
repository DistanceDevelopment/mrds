#' beta
#' 
#' this function derives the integral of the density divided by key and adj at
#' distance zero
#' 
#' look at int1 <- integratepdf(ddfobj, select=!x$binned, width=right,
#' int.range=intrange[!x$binned, ],
#' doeachint=doeachint,
#' point=misc.options$point, standardize=standardize,
#' left=left)
#' in flpt.lnl.r for inspiration on how to use this function. 

beta <- function(int.range, ddfobj, width, standardize = TRUE, 
                 point = FALSE, left = 0, stdint = FALSE, select = NULL,
                 index = 1) {
  
  ## derive the integral of g
  g <- gstdint(int.range[1, ], ddfobj = ddfobj, index = index, select = select,
               width = width, standardize = standardize, point = point,
               stdint = stdint, left = left)
  
  ## Derive adj term and series, and key for distance = 0
  key.val.0 <- switch(key,
                      hn    = keyfct.hn(0, key.scale),
                      hr    = keyfct.hz(0, key.scale, key.shape),
                      unif  = 1,
                      gamma = keyfct.gamma(0, key.scale, key.shape),
                      th1   = keyfct.th1(0, key.scale, key.shape),
                      th2   = keyfct.th2(0, key.scale, key.shape))
  
  adj.val.0 <- switch(adj.series,
                      poly = adjfct.poly(0, scaling, adj.order,
                                         adj.parm, adj.exp),
                      herm = adjfct.herm(0, scaling, adj.order,
                                         adj.parm, adj.exp),
                      cos  = adjfct.cos(0, scaling, adj.order,
                                        adj.parm, adj.exp))
  
  return(g / (key.val.0 * (1 + adj.val.0)))
}