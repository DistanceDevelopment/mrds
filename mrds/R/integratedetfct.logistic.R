integratedetfct.logistic <-
function(x,scalemodel,width,theta1,integral.numeric,w)
#
# integratedetfct.logistic - integrates a logistic detection function; a separate function is used because
# in certain cases the integral can be solved analytically and also because the scale trick used with the half-normal
# and hazard rate doesn't work with the logistic.
#
# Arguments:
#
# x                - logistic design matrix values
# scalemodel       - scale model for logistic
# width            - transect width
# theta1           - parameters for logistic
# integral.numeric - if TRUE computes numerical integral value
# w                - design covariates
#
# Value:
# 
# vector of integral values
# 
# Functions Used:  is.logistic.constant, integratelogistic, integratelogistic.analytic
# 
{
#
#     If the integrals are to be computed numerically use integratelogistic function
#
      if(integral.numeric)
      {
#
#         If all values are constant (no different covariates) then repeat integral value for a single
#         integration
#
          if(is.logistic.constant(x,scalemodel,width))
          {
             int1 <- rep(integratelogistic (x=x[1,], models=list(g0model=scalemodel), beta=theta1,width),dim(x)[1]) 
          } else
#
#         Otherwise compute integrals for all different sets of covariates
#
          {
             int1 <- NULL
             for (i in 1:(dim(x)[1]))
                  int1 <- c(int1, integratelogistic (x=x[i,], models=list(g0model=scalemodel), beta=theta1,width))
          }
#
#     If not a numerical integral use integratelogistic.analytic to compute analytically
#
      } else
          int1 <- integratelogistic.analytic(x, models=list(g0model=scalemodel), beta=theta1, width=width)
#
#     Scale of logistic at x=0 to scale detection function to g(0)=1
#
      int1=int1/g0(beta=theta1,w[[2]])
      return(int1)
}
