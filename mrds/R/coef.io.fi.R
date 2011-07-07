
coef.io.fi <-
function(object,...)
{
#
# coef.io.fi
#
# Provides a summary of parameters and estimates from the output of io.fi object
#
# Arguments:
#
# object      - object from ddf.io.fi
#
# Value: list of coefficient data frames (scale and exponent (if hazard))
#
   if(length(grep("gam",as.character(object$model)))==0)
      vcov=solvecov(object$hessian)$inv
   else
      vcov=object$hessian
   coeff=as.data.frame(cbind(estimate=coef(object$mr),se=sqrt(diag(vcov))))
   return(coeff)
}
