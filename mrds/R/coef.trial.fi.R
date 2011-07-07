coef.trial.fi <-
function(object,...)
{
#
# coef.trial.fi
#
# Provides a summary of parameters and estimates from the output of trial.fi object
#
# Arguments:
#
# object      - object from ddf.trial.fi
#
# Value: list of coefficient data frames (scale and exponent (if hazard))
#
   vcov=solvecov(object$hessian)$inv
   coeff=as.data.frame(cbind(estimate=coef(object$mr),se=sqrt(diag(vcov))))
   return(coeff)
}
