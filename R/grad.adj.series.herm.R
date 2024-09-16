#' Series of the gradient of the Hermite polynomial adjustment series w.r.t. the 
#' scaled distance.
#' 
#' For internal use only -- not to be called by 'mrds' or 'Distance' users 
#' directly.
#'
#' @param distance perpendicular distance vector/scalar
#' @param scaling scale parameter
#' @param adj.order the adjustment order
#' @param adj.parm vector of parameters (a_j)
#' @param adj.exp boolean, defaults to FALSE
#' 
#' @return scalar or vector containing the gradient of the Hermite adjustment 
#' series for every value in \code{distance} argument
#'
#' @author Felix Petersma
grad.adj.series.herm <- function(distance, scaling = 1, adj.order, 
                                 adj.parm = NULL, adj.exp = FALSE){
  
  # Check the adjustment parameters
  if(is.null(adj.parm)){
    adj.parm <- as.vector(rep(1, length(adj.order)))
  }
  
  adj.order <- as.vector(adj.order)
  
  hermsum <- 0
  
  for(i in seq_along(adj.order)){
    if (adj.order[i] - 1 == 0) {
      hermsum <- hermsum + 
        (adj.parm[i] * adj.order[i] * 1)
    } else {
      hermsum <- hermsum + 
        (adj.parm[i] * adj.order[i] * hermite.poly((distance / scaling),
                                                   adj.order[i] - 1))
    }
  }
  
  # if adj.exp return exp(hermsum) to keep f(x)>0
  if(adj.exp){
    return(exp(hermsum))
  }else{
    return(hermsum)
  }
}
