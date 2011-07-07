logisticbyx <-
function (distance, x, models, beta, point)
#
#  logisticbyx - treats logistic as a function of covariate z; for a given z it computes
#                function at with those covariate values at a range of distances 
#
#  Arguments:
#  
#  distance - vector of distance values
#  x        - covariate data
#  models   - model list
# 
#  value: vector of probabilities
#  
#  Functions used: g0, setcov
# 
{
   xlist <- as.list(x)
   xlist$distance <-distance
   xmat <- expand.grid(xlist)
   if(!point)
      return(g0(beta, setcov(xmat, models$g0model)$cov))
   else
	   return(g0(beta, setcov(xmat, models$g0model)$cov)*2*distance)
   
}
