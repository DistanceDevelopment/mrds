logisticdetfct <-
function(distance, theta, w, std=FALSE)
#
# logisticdetfct
#
# Logistic detection function with x = distance and z = other covariates, X(x,z) is 
# design matrix for x,z.  
#
# h(x,z) = exp (X(x,z)*beta)/(1+exp(X(x,z)*beta)
# g(x,z)= h(x,z)/h(0,z)
#
# Arguments:
#
# distance - perpendicular distance vector
# theta    - scale parameters
# w        - scale covariate matrix 
# std      - if TRUE uses scale=1
#
# Value:
#
# The routine returns a vector of probabilities that the observation
# were detected given they were at the specified distance and assuming that
# g(0)=1 (ie a standard line transect detection function).
#
{
  if(is.matrix(w[[1]])) 
     exp(w[[1]] %*% theta)/(1 + exp(w[[1]] %*% theta))/g0(beta=theta,w[[2]])
  else
     exp(as.matrix(w[[1]]) %*% beta)/(1 + exp(as.matrix(w[[1]]) %*% beta))/g0(beta=theta,w[[2]])
}
