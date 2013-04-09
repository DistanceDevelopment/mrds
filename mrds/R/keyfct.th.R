keyfct.th1 <- function(distance, key.scale, key.shape)
{
#
# keyfct.th1
#
# Threshold key function
#
# Arguments:
# distance  - perpendicular distance vector
# key.scale - vector of scale values
# key.shape - vector of shape values
#
# Value: vector of probabilities  
#
	
	return( 0.5 - 0.5*erf(distance/key.scale - key.shape))
}
keyfct.th2 <- function(distance, key.scale, key.shape)
{
#
# keyfct.th2
#
# Threshold key function
#
# Arguments:
# distance  - perpendicular distance vector
# key.scale - vector of scale values
# key.shape - vector of shape values
#
# Value: vector of probabilities  
#
	return(erf(exp(key.shape-distance/key.scale)))
}
# error function
erf <- function(x) 2*pnorm(x*sqrt(2)) - 1
