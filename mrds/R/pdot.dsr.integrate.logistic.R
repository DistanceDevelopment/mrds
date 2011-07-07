pdot.dsr.integrate.logistic <-
function(right, width, beta, x, integral.numeric, BT, models, GAM=FALSE, rem=FALSE, point=FALSE) 
{
# 
#   pdot.dsr.integrate.logistic
#
#   Compute probability that it was detected from at least one observer (pdot)
#   for a logistic g* that contains distance
# 
#   Functions called:
#
#   integratelogistic- computes integral (numerically) over x from 0 to width of a 
#                        logistic detection function 
#   integratelogisticdup- computes integral (numerically) over x from 0 to width of a 
#                          duplicate logistic detection function 
#   integratelogistic.analytic- computes integral (analytically) over x from 0 to width of a 
#                          logistic detection function 
#   integrate  - S+ routine to integrate a function 
#   logisticbyx- computes the detection probability with a logistic det fct for a single set of 
#                                   z covariates but multiple values of x
#   logisticdupbyx- computes the duplicate detection probability with a logistic det fct 
#                       for a single set of z covariates but multiple values of x
#   logisticbyz- computes the detection probability at a single x with a logistic det fct 
#             with mulitple sets of z covariates 
#
#
#   Uniform detection function for g' but g* includes distance   
#
#
#  If the models are non-linear in distance, numerical integation is required for int1 and int2
#  
   if(integral.numeric | point)
   {
      if(GAM)
         is.constant=FALSE
      else
         is.constant=is.logistic.constant(x[x$observer==1,],models$g0model,width)
      if(is.constant)
      {
         int1 <- rep(integratelogistic (x=(x[x$observer==1,])[1,], models, beta,right, point),dim(x[x$observer==1,])[1]) 
      } else
      {
        int1 <- NULL
        for (i in 1:(dim(x[x$observer==1,])[1]))
           int1 <- c(int1, integratelogistic (x=(x[x$observer==1,])[i,], models, beta,right,point))
      }
      if(!BT)
      {
         if(is.logistic.constant(x[x$observer==2,],models$g0model,width))
         {
            int2 <- rep(integratelogistic (x=(x[x$observer==2,])[1,], models, beta,right, point),
                              dim(x[x$observer==2,])[1]) 
         } else
         {
            int2 <- NULL
            for (i in 1:(dim(x[x$observer==2,])[1]))
                 int2 <- c(int2, integratelogistic (x=(x[x$observer==2,])[i,], models, beta,right, point))
         }
      } 
      else
          int2 <- NULL
   } else
#
#  If the models are linear in distance, solve int1 and int2 analytically 
#
   {
      int1 <- integratelogistic.analytic(x[x$observer==1,], models=models, beta=beta, width=right)
      if(!BT) 
          int2 <- integratelogistic.analytic(x[x$observer==2,], models=models, beta=beta, width=right)
      else
           int2 <- NULL
   }
#
#  Numerical integration is always required for int3
#
   if(!BT )
   {
     if(is.logistic.constant(x[x$observer==1,],models$g0model,width) & 
        is.logistic.constant(x[x$observer==2,],models$g0model,width))
     {
         int3 <- rep(integratelogisticdup(x1=(x[x$observer==1,])[1,],
                 x2=(x[x$observer==2,])[1,],models,beta,right, point),dim(x[x$observer==2,])[1])
     } else
     {
         int3 <- NULL
         for (i in 1:(dim(x[x$observer==1,])[1]))
             int3 <- c(int3, integratelogisticdup(x1=(x[x$observer==1,])[i,],x2=(x[x$observer==2,])[i,],models,beta,right, point))
     }
     pdot <- int1 + int2 - int3
   } else
   {
     int3 <- NULL
     pdot <- int1
   }

   if(!point)
      div=width
  else 
      div=width^2
#
#   Return list of results
#
   return(list(pdot=pdot/div,int1=int1/div,int2=int2/div,int3=int3/div))
}
