#' Fit detection function using key-adjustment functions
#'
#' Fit detection function to observed distances
#' using the key-adjustment
#' function approach.  If adjustment functions are included it will alternate
#' between fitting parameters of key and adjustment functions and then all
#' parameters much like the approach in the CDS and MCDS Distance FORTRAN code.
#' This function is called by the driver functioin \code{detfct.fit}.  This 
#' function does the calls the optimx() function (from the package optimx).
#'
#' @import optimx Rsolnp
#' @aliases detfct.fit.opt
#' @param ddfobj detection function object
#' @param optim.options control options for optim
#' @param bounds bounds for the parameters
#' @param misc.options miscellaneous options
#' @param fitting character string with values "all","key","adjust" to
#'   determine which parameters are allowed to vary in the fitting
#' @return fitted detection function model object with the following list
#'   structure \item{par}{final parameter vector} \item{value}{final negative
#'   log likelihood value} \item{counts}{number of function evaluations}
#'   \item{convergence}{see codes in optim} \item{message}{string about
#'   convergence} \item{hessian}{hessian evaluated at final parameter values}
#'   \item{aux}{ a list with 20 elements \itemize{ \item maxit: maximum number
#'   of iterations allowed for optimization \item lower: lower bound values for
#'   parameters \item upper: upper bound values for parameters \item setlower:
#'   TRUE if they are user set bounds \item setupper: TRUE if they are user set
#'   bounds \item point: TRUE if point counts and FALSE if line transect \item
#'   int.range: integration range values \item showit: integer value that
#'   determines information printed during iteration \item integral.numeric
#'   if TRUE compute logistic integrals numerically \item
#'   breaks: breaks in distance for defined fixed bins for analysis \item
#'   maxiter: maximum iterations used \item refit: if TRUE, detection function
#'   will be fitted more than once if parameters are at a boundary or when
#'   convergence is not achieved \item nrefits: number of refittings \item
#'   parscale: parameter scale values
#'   \item mono: if TRUE, montonicity will be enforced \item
#'   mono.strict: if TRUE, then strict monotonicity is enforced; otherwise weak
#'   \item width: radius of point count or half-width of strip \item
#'   standardize: if TRUE, detection function is scaled so g(0)=1 \item ddfobj:
#'   distance detection function object; see \code{\link{create.ddfobj}} \item
#'   bounded: TRUE if parameters ended up a boundary (I think) \item model:
#'   list of formulas for detection function model (probably can remove this)
#'   }}
#' @author Dave Miller; Jeff Laake; Lorenzo Milazzo
detfct.fit.opt <- function(ddfobj, optim.options, bounds, misc.options,
                           fitting="all"){
  # Functions Used: assign.par, detfct.fit.opt, errors, get.par

  # grab the initial values
  initialvalues <- getpar(ddfobj)
  initialvalues.set <- initialvalues # store for later

  bounded <- TRUE
  refit.count <- 0
  # Set some shortcuts
  lowerbounds <- bounds$lower
  upperbounds <- bounds$upper
  showit <- misc.options$showit
  refit <- misc.options$refit
  nrefits <- misc.options$nrefits

  # jll 18-sept-2006; added this code to get the logicals that indicate whether
  # lower/upper bound settings were specified by the user
  setlower <- bounds$setlower
  setupper <- bounds$setupper

  # grab the method(s) if we're using optimx()
  if(!misc.options$mono){
    opt.method <- optim.options$optimx.method
    optim.options$optimx.method <- NULL
    optim.options$follow.on <- TRUE
  }

  # if monotonicity has been requested but we are using key only then just
  # use optimx
  if(misc.options$mono & is.null(ddfobj$adjustment)){
    misc.options$mono <- FALSE
  }

  # if nofit=TRUE, we just want to set the parameters, calculate the
  # likelihood and exit
  if(misc.options$nofit){
    if(!is.null(initialvalues) && any(is.na(initialvalues))){
      stop("No fitting, but initial values not specified!\n")
    }

    lt <- list()
    lt$par <- initialvalues
    lt$value <- flnl(initialvalues,ddfobj,misc.options)
    lt$hessian <- NULL
    lt$model <- list(scalemodel=misc.options$scalemodel)
    lt$converge <- 0
    lt$message<-"MAYBE CONVERGENCE?"

    lt$aux <- c(optim.options,bounds,misc.options)
    ddfobj <- assign.par(ddfobj,lt$par)
    lt$aux$ddfobj <- ddfobj

    lt$optim.history<-rbind(misc.options$optim.history,
                            c(lt$conv,-lt$value,lt$par))

    lt$conv<-NULL
    return(lt)
  }
  # 30 Jan 06; jll- modified default parscale values
  if(any(is.na(misc.options$parscale))){
    misc.options$parscale <- abs(initialvalues)
  }else{
    if(length(misc.options$parscale)!=length(initialvalues)){
      errors("Incorrect length of parameter scale vector; using default values\n")
      misc.options$parscale <- abs(initialvalues)
    }
  }

  # save last value of the lnl -- starting value
  lnl.last <- Inf

  # recover the optimisation history
  optim.history <- misc.options$optim.history

  # while parameters are bounded continue refitting and adjusting bounds
  while(bounded){
    itconverged <- FALSE

    # Continue fitting until convergence occurs as long as refitting
    # is requested
    while(!itconverged){
      # Call optimization routine to find constrained mles; upon
      # completion add the user specified models and return the list.

      # if we want monotonicity, use Lorenzo's code...
      if(misc.options$mono){
        # lower and upper bounds of the inequality constraints
        lowerbounds.ic <- rep(0,2*misc.options$mono.points)
        upperbounds.ic <- rep(10^10,2*misc.options$mono.points)

        #lt<-try(solnp(pars=initialvalues, fun=flnl, eqfun=NULL, eqB=NULL,
        #              ineqfun=flnl.constr,
        #              ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
        #              LB=lowerbounds, UB=upperbounds,
        #              ddfobj=ddfobj, misc.options=misc.options,
        #              control=list(trace=as.integer(showit),
        #                           tol=misc.options$mono.tol,
        #                           delta=misc.options$mono.delta)))

        # this code randomly generates starting values see ?gosolnp
        lt<-try(gosolnp(pars=initialvalues, fun=flnl, eqfun=NULL, eqB=NULL,
                      ineqfun=flnl.constr,
                      ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                      LB=lowerbounds, UB=upperbounds,
                      ddfobj=ddfobj, misc.options=misc.options,
                      control=list(trace=as.integer(showit),
                                   tol=misc.options$mono.tol,
                                   delta=misc.options$mono.delta),
                      distr = rep(1, length(lowerbounds)),
                      n.restarts = 2, n.sim = 200,
                      rseed=as.integer(runif(1)*1e9)))


        # if that failed then make a dummy object
        if(class(lt)=="try-error"){
          lt <- list()
          lt$conv <- 9
          lt$value <- lnl.last
          lt$par <- initialvalues

          if(showit==3){
            errors("Optimisation failed, ignoring and carrying on...")
          }
        }else{
          # above gosolnp code stores best lnl as last value
          # recover that information!
          lt$conv <- lt$convergence[length(lt$values)]
          lt$par <- lt$pars[length(lt$values)]
          lt$value <- lt$values[length(lt$values)]
          lt$message <- ""
        }

        # re-jig some stuff so that this looks like an optim result...
        lt$conv <- lt$convergence
        lt$par <- lt$pars
        lt$value <- lt$values[length(lt$values)]
        lt$message <- ""

      }else{
        # use optimx
        lt <- try(optimx(initialvalues, flnl, method=opt.method,
                     control=optim.options,
                     hessian=TRUE, lower=lowerbounds,
                     upper=upperbounds,ddfobj=ddfobj, fitting=fitting,
                     misc.options=misc.options))

        topfit.par <- coef(lt, order="value")[1, ]
        details <- attr(lt,"details")[1,]
        lt <- as.list(summary(lt, order="value")[1, ])
        lt$par <- topfit.par
        lt$message <- ""
        names(lt)[names(lt)=="convcode"] <- "conv"
        lt$hessian <- details$nhatend
      }

      # Print debug information
      if(showit==3){
        errors(paste("Converge = ",lt$conv,"\n",
                     "lnl = ",lt$value,"\n",
                     "parameters = ",paste(lt$par,collapse=", ")))
      }

      optim.history<-rbind(optim.history,c(lt$conv,-lt$value,lt$par))

      ## Convergence?!
      # OR... don't wiggle the pars or do refits if we're doing adjustment
      # or key fitting alone only do it in "all" mode
      if(lt$conv==0|!refit | fitting!="all"){
        itconverged<-TRUE

        lt$aux <- c(optim.options,bounds,misc.options)
        ddfobj <- assign.par(ddfobj,lt$par)
        lt$aux$ddfobj <- ddfobj
      }else{
      # If we don't have convergence what do we do
        refit.count<-refit.count+1
        if(is.null(nrefits)|refit.count<=nrefits){
          if(showit>=1){
            errors("No convergence. Refitting ...")
          }

          # if the new values weren't as good, take the last set
          # and jiggle them a bit...

          # previously this was a bit weird? That addition should be the
          # same sign (and be allowed to be negative) otherwise we're
          # just making it more positive
          #initialvalues <- lt$par*(runif(length(initialvalues))+.5)

          initialvalues <- lt$par*runif(length(initialvalues),
                                        sign(lowerbounds-1),
                                        sign(upperbounds+1))

          lnl.last <- lt$value

          # if we just replace NAs with 0s then sometimes these are out
          # of bounds so replace with the values we started with
          initialvalues[is.na(initialvalues)] <-
                                      initialvalues.set[is.na(initialvalues)]
        }else{
          itconverged <- TRUE
        }
      }
    }

    if(any(is.na(lt$par)) | lt$conv!=0){
      # if there was no convergence then just return the lt object for debugging
      errors("Problems with fitting model. Did not converge")
      if(misc.options$debug){
        lt$optim.history <- optim.history
        return(lt)
      }else{
        stop("No convergence.")
      }
    }

    # check whether parameters hit their bounds
    bounded <- check.bounds(lt,lowerbounds,upperbounds,ddfobj,
                            showit,setlower,setupper)

    if(!refit){
        bounded <- FALSE
    }else{
      refit.count <- refit.count+1
      if(!is.null(nrefits)){
        if(refit.count>nrefits){
          bounded <- FALSE
          lt$message <- "FALSE CONVERGENCE"
          lt$conv <- 1
        }
      }
    }

    # fix: jll 18-Nov-04; previous code would get stuck at 0 if sign of
    # initial value was opposite of the mle.
    # Additional statement skips over 0.
    if(bounded){
      bound.low <- abs(lt$par-lowerbounds)<1e-6
      bound.hi <- abs(lt$par-upperbounds)<1e-6
      if(!setlower) {
        lowerbounds[bound.low] <- lowerbounds[bound.low] -
                                   0.5*abs(lowerbounds[bound.low])
        lowerbounds[bound.low & lowerbounds>0 &lowerbounds < 0.5] <- -0.5
      }
      if(!setupper){
        upperbounds[bound.hi] <- upperbounds[bound.hi] +
                                  0.5*abs(upperbounds[bound.hi])
        upperbounds[bound.hi & upperbounds<0 &upperbounds > -0.5] <- 0.5
      }

      if(showit>=1){
        errors("Refitting ...")
      }

    }
  }
  lt$model <- list(scalemodel=misc.options$scalemodel)
  lt$converge <- lt$conv
  lt$conv <- NULL

  # save bounded status
  lt$bounded <- bounded

  # save the bounds
  bounds$lower <- lowerbounds
  bounds$upper <- upperbounds
  lt$bounds <- bounds

  # save optimisation history
  lt$optim.history <- optim.history

  return(lt)
}
