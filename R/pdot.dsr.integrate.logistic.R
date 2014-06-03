#' Compute probability that a object was detected by at least one observer
#'
#' Computes probability that a object was detected by at least one observer
#' (\code{pdot} or p_.) for a logistic detection function that contains
#' distance.
#'
#' @param right ???
#' @param width transect width
#' @param beta parameters of logistic detection function
#' @param x data matrix
#' @param integral.numeric ???
#' @param BT ???
#' @param models list of models including \code{g0model}
#' @param GAM ???=FALSE
#' @param rem ???=FALSE
#' @param point \code{TRUE} for point transects
#'
#' @author Jeff Laake
#'
pdot.dsr.integrate.logistic <- function(right, width, beta, x,
                                        integral.numeric, BT, models, GAM=FALSE,
                                        rem=FALSE, point=FALSE){
  # Functions called:
  # integratelogistic - computes integral (numerically) over x from 0 to width
  #                     of a logistic detection function
  # integratelogisticdup - computes integral (numerically) over x from 0 to
  #                        width of a duplicate logistic detection function
  # integratelogistic.analytic - computes integral (analytically) over x from 0
  #                              to width of a logistic detection function
  # logisticbyx - computes the detection probability with a logistic det fct
  #               for a single set of z covariates but multiple values of x
  # logisticdupbyx - computes the duplicate detection probability with a
  #                  logistic det fct for a single set of z covariates but
  #                  multiple values of x
  # logisticbyz - computes the detection probability at a single x with a
  #               logistic det fct with mulitple sets of z covariates
  #
  # Uniform detection function for g' but g* includes distance
  #
  # If the models are non-linear in distance, numerical integation is
  # required for int1 and int2

  if(length(right)>1){
    lower <- right[1]
    right <- right[2]
  }else{
    lower <- 0
  }

  if(is.null(x$observer)){
    stop("data must have a column named \"observer\"")
  }

  if(integral.numeric | point){

    if(GAM|length(right)>1){
      is.constant <- FALSE
    }else{
      is.constant <- is.logistic.constant(x[x$observer==1,],
                                          models$g0model,width)
    }

    if(is.constant){
      int1 <- rep(integratelogistic(x=x[x$observer==1,][1,], models, beta,
                                    lower=0,right, point),
                  nrow(x[x$observer==1,]))
    }else{
      int1 <- NULL
      for(i in 1:nrow(x[x$observer==1,])){
        int1 <- c(int1, integratelogistic(x=(x[x$observer==1,])[i,], models,
                                          beta,lower=lower,right,point))
      }
    }

    if(!BT){
      if(is.logistic.constant(x[x$observer==2,],models$g0model,width)){
        int2 <- rep(integratelogistic(x=x[x$observer==2,][1,], models, beta,
                                      lower=lower,right, point),
                    nrow(x[x$observer==2,]))
      }else{
        int2 <- NULL
        for(i in 1:nrow(x[x$observer==2,])){
          int2 <- c(int2, integratelogistic(x=x[x$observer==2,][i,], models,
                                            beta,lower=lower,right, point))
        }
      }
    }else{
      int2 <- NULL
    }
  }else{
    # If the models are linear in distance, solve int1 and int2 analytically
    int1 <- integratelogistic.analytic(x[x$observer==1,], models=models,
                                       beta=beta, width=right)
    if(!BT){
      int2 <- integratelogistic.analytic(x[x$observer==2,], models=models,
                                         beta=beta, width=right)
    }else{
      int2 <- NULL
    }
  }

  # Numerical integration is always required for int3
  if(!BT){
    if(is.logistic.constant(x[x$observer==1,],models$g0model,width) &
       is.logistic.constant(x[x$observer==2,],models$g0model,width)){

      int3 <- rep(integratelogisticdup(x1=x[x$observer==1,][1,],
                                       x2=x[x$observer==2,][1,],models,beta,
                                       lower=lower,right, point),
                  nrow(x[x$observer==2,]))
    }else{
      int3 <- NULL
      for(i in 1:nrow(x[x$observer==1,])){
        int3 <- c(int3, integratelogisticdup(x1=(x[x$observer==1,])[i,],
                                             x2=(x[x$observer==2,])[i,],
                                             models, beta, lower=lower,
                                             right, point))
      }
    }
    pdot <- int1 + int2 - int3
  }else{
    int3 <- NULL
    pdot <- int1
  }

  if(!point){
    div <- width
  }else{
    div <- width^2
  }

  # Return list of results
  return(list(pdot = pdot/div,
              int1 = int1/div,
              int2 = int2/div,
              int3 = int3/div))
}
