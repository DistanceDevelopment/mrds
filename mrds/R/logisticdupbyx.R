logisticdupbyx <-
function (distance, x1, x2, models, beta, point)
#
#  logisticdupbyx - treats logistic for duplicates as a function of covariate z; for a given z it computes
#                   function at with those covariate values at a range of distances 
#
#  Arguments:
#  
#  distance - vector of distance values
#  x1       - covariate data for fct 1
#  x2       - covariate data for fct 2
#  models   - model list
#  beta     - parameters 
#
#  value: vector of probabilities
#  
#  Functions used: g0, setcov
{
   xlist <- as.list(x1)
   xlist$distance <-distance
   xmat <- expand.grid(xlist)
   gx1 <- g0(beta, setcov(xmat, models$g0model)$cov)
   xlist <- as.list(x2)
   xlist$distance <-distance
   xmat <- expand.grid(xlist)
   if(!point)
	   return(gx1 * g0(beta, setcov(xmat, models$g0model)$cov))
   else
	   return(gx1 * g0(beta, setcov(xmat, models$g0model)$cov)*2*distance)   
}
