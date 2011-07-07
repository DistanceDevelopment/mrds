g0 <-
function(beta, z)
#
#  g0 - computes value of p(0) using a logit formulation
#
#  Arguments:
#
#  beta - parameters
#  z    - design matrix of covariate values
# 
#  value: vector of p(0) values
#
#
{
#
#  Treat differently if z is an intercept vector or a design matrix
#
  if(is.matrix(z)) 
     exp(z %*% beta)/(1 + exp(z %*% beta))
  else
     exp(as.matrix(z) %*% beta)/(1 + exp(as.matrix(z) %*% beta))    
}
