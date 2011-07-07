is.logistic.constant <-
function(xmat,g0model,width)
#
#  is.logistic.constant - determines whether the specified logit model is constant for all observations 
#  If it is constant then only one integral needs to be computed
#
# Arguments:
#
# xmat     - data
# g0model  - logit model
# width    - transect width
#
# Value: logical value 
#
{ 
  xmat$distance <- rep(width, dim(xmat)[1])
    zlist <-setcov(xmat, g0model)
    beta <- rep(1,zlist$dim)
  logit1 <- beta %*% t(zlist$cov)
    if(all(logit1[1]==logit1))
       return(TRUE)
    else
       return(FALSE)
}
