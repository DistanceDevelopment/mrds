adjfct.cos <-
		function(distance,scaling=1,adj.order,adj.parm=NULL)
#
# adjfct.cos
#
# Cosine adjustment terms
#
#
# Arguments:
#
# distance  - perpendicular distance vector
# scaling     - scale parameter
# adj.order - vector of orders of Cosine terms to fit
# adj.parm  - vector of parameters (a_j)
#
{
# Check the adjustment parameters
	if(is.null(adj.parm)){
		adj.parm <- as.vector(rep(1,length(adj.order)))
	}
	
	adj.order <- as.vector(adj.order)
	
	cossum <- 0
	
	for(i in 1:length(adj.order))
		cossum <- cossum + (adj.parm[i]*cos((adj.order[i]*pi*distance)/scaling))
	
	### EXPERIMENTAL "feature"
# To try and keep the likelihood returning a postive value (so that it may be log'd)
# we make the cossum > -1 
	
#  if(cossum == -1)
#    cossum <- 0
	
	return(cossum)
}
