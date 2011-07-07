scalevalue <-
function(key.scale, z)
#
# scalevalue - computes scale of detection function for a given set of scale covariates and parameters
#              It uses a log link
# 
# Arguments:
#
#  key.scale - scale parameters
#  z	     - design matrix for scale covariates
#
# Value:
#     
# Vector of scale values
#
{
  if(is.matrix(z)) 
    exp(z %*% key.scale)
  else
    exp(as.matrix(z) %*% key.scale)
}
