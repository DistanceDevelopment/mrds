#' Generate data from a fitted detection function
#'
#' @param ds.object a fitted detection function object
#'
#' @note This function changes the random number generator seed. To avoid any potential side-effects, use something like: \code{seed <- get(".Random.seed",envir=.GlobalEnv)} before running code and \code{assign(".Random.seed",seed,envir=.GlobalEnv)} after.
#' @author David L. Miller
#' @importFrom stats runif rnorm fitted
sample_ddf <- function(ds.object){

  # this would make life a lot easier...
  if(!inherits(ds.object, "ds")){
    stop("Can only sample from ds models")
  }

  # how many samples do we need?
  n.ds.samples <- nrow(ds.object$data)

  # how many samples do we have so far?
  n.samps <- 0
  dists <- c()

  # make sure that a model gets fitted
  dud.df <- TRUE

  covars <- FALSE
  direct_sim <- FALSE
  # can we directly simulate?
  if(ds.object$ds$aux$ddfobj$type=="hn" &
     is.null(ds.object$ds$aux$ddfobj$adjustment$parameters) &
     ds.object$ds$aux$ddfobj$scale$formula=="~1"){
    direct_sim <- TRUE
  }
  # do we have a covariate model?
  if(ds.object$ds$aux$ddfobj$scale$formula != "~1"){
    covars <- TRUE
    formula <- as.formula(ds.object$ds$aux$ddfobj$scale$formula)
    model_xmat <- ds.object$ds$aux$ddfobj$xmat
  }

  # create an object to hold the parameters
  pars <- list()
  pars$scale <- ds.object$ds$aux$ddfobj$scale$parameters
  if(!is.null(ds.object$ds$aux$ddfobj$shape$parameters)){
    pars$shape <- ds.object$ds$aux$ddfobj$shape$parameters
  }
  if(!is.null(ds.object$ds$aux$ddfobj$adjustment$parameters)){
    pars$adjustment <- ds.object$ds$aux$ddfobj$adjustment$parameters
  }


  # okay now do the simulation
  while(dud.df){
    if(direct_sim){
      # if we just have a half-normal key function then we can
      # directly simulate...
      dists <- abs(rnorm(n.ds.samples, mean=0, sd=exp(pars$scale)))
      dists <- data.frame(distance = dists,
                          detected = rep(1,length(dists)),
                          object   = 1:length(dists))

    }else{
      # otherwise we need to do some rejection sampling

      # what should the sampler be?
      if(ds.object$meta.data$point){
        sampler <- function(n) rtriangle(n, a=ds.object$meta.data$left,
                                         b=ds.object$meta.data$width,
                                         c=ds.object$meta.data$width)
      }else{
        sampler <- function(n) runif(n, min=ds.object$meta.data$left,
                                     max=ds.object$meta.data$width)
      }
      # what function do we use to get the acceptance probability
      accept_p <- mrds::detfct

      # since rejection sampling is time consuming, generate lots of
      # samples at once, we re-scale the number by the inverse of the 
      # ratio accepted. The first time over, let's make that 10x
      n.mult <- 10

      # data
      xdat <- ds.object$data
      xdat$distance <- rep(NA, nrow(xdat))
      xdat$binned <- rep(FALSE, nrow(xdat))
      # condition on covariates
      if(covars){
        vars <- all.vars(formula)
        xdat[, vars] <- model_xmat[, vars]
      }

      while(n.samps < n.ds.samples){

        # which samples should we take?
        #if(is.null(dists$object)){
          ind <- xdat$object
        #}else{
        #  ind <- xdat$object[!(xdat$object %in% dists$object)]
        #}
        this.n.samps <- n.mult*length(ind)

        # generate some new distances
        this_xdat <- xdat[rep(which(ind %in% xdat$object), n.mult),]
        this_xdat$distance <- sampler(this.n.samps)


        # create a ddf object
        ddfobj <- mrds::create.ddfobj(as.formula(ds.object$dsmodel), this_xdat,
                                      ds.object$meta.data, pars)

        # generate acceptance probability
        U <- runif(this.n.samps)

        # was it accepted?
        inout <- U <= accept_p(this_xdat$distance, ddfobj,
                               standardize=FALSE, width=ds.object$meta.data$width)

        dists <- rbind(dists, this_xdat[inout, ])
        dists <- dists[!duplicated(dists$object), ]

        n.samps <- nrow(dists)

        # update the number of extra samples we make by inverting the ratio
        # of accepted to generated this round
        n.mult <- min(1, ceiling(1/(sum(inout)/this.n.samps)))
      }
    }

    # make sure that we got the right number
    dists <- dists[!duplicated(dists$object), ]

    # fit the model to the new data
    ddf.call <- ds.object$call
    ddf.call$data <- dists
    ddf.call$meta.data <- ds.object$meta.data
    ddf.call$control <- ds.object$control
    ddf.call$control$initial <- pars
    ddf.call$dsmodel <- as.formula(ds.object$dsmodel)
    ddf.fitted <- suppressWarnings(suppressMessages(try(
                    with(ds.object, eval(ddf.call)))))

    # if it all went well, then set dud.df to FALSE and quit the loop
    if(all(class(ddf.fitted) != "try-error")){
      dud.df <- FALSE
    }else{
      # otherwise forget everything and start again
      n.samps <- 0
      dists <- c()
    }
  }

  # return the offset
  return(ddf.fitted)
}
