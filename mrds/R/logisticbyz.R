logisticbyz <-
function (x, distance, models, beta)
#
#  logisticbyz - treats logistic as a function of distance x; for a given x it computes
#                function at all covariate values in data
#
#  Arguments:
#
#  x        - covariate data  
#  distance - single distance value
#  models   - model list
#  beta     - model parameters 
#
#  value: vector of probabilities
#  
#  Functions used: g0, setcov
# 

{
    x$distance <- rep(distance,length(x$distance))
    zlist <- setcov(x, models$g0model)
    g0(beta,zlist$cov)
}
