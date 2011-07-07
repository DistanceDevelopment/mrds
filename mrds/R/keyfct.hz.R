keyfct.hz <-
		function(distance, key.scale, key.shape)
#
# keyfct.hz
#
# Hazard rate key function: 1 - exp (- (x/scale) ^ (-shape)))
#
# Arguments:
# distance  - perpendicular distance vector
# key.scale - vector of scale values
# key.shape - vector of shape values
#
# Value: vector of probabilities that the observations were detected 
# given they were at the specified distance and assuming that g(0)=1 
# (ie a standard line transect detection function).
#
{
	return(
			1 - exp( - (distance/key.scale)^( - key.shape))
	)
}

