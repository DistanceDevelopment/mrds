#' Fit detection function using key-adjustment functions
#'
#' Fit detection function to observed distances using the key-adjustment
#' function approach. If adjustment functions are included it will alternate
#' between fitting parameters of key and adjustment functions and then all
#' parameters much like the approach in the CDS and MCDS Distance FORTRAN code.
#' This function is called by the driver function \code{detfct.fit}, it then
#' calls the relevant optimisation routine, \code{\link[nloptr]{slsqp}},
#' \code{\link[Rsolnp]{solnp}} or \code{\link[optimx]{optimx}}.
#'
#' @import nloptr optimx Rsolnp
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
#'   convergence is not achieved \item nrefits: number of refittings
#'   \item mono: if TRUE, monotonicity will be enforced \item
#'   mono.strict: if TRUE, then strict monotonicity is enforced; otherwise weak
#'   \item width: radius of point count or half-width of strip \item
#'   standardize: if TRUE, detection function is scaled so g(0)=1 \item ddfobj:
#'   distance detection function object; see \code{\link{create.ddfobj}} \item
#'   bounded: TRUE if estimated parameters are at the bounds \item model:
#'   list of formulas for detection function model (probably can remove this)
#'   }}
#' @author Dave Miller; Jeff Laake; Lorenzo Milazzo; Felix Petersma
#' @importFrom stats runif optim
detfct.fit.opt <- function(ddfobj, optim.options, bounds, misc.options,
                           fitting="all"){
  
  # grab the initial values
  initialvalues <- getpar(ddfobj)
  initialvalues.set <- initialvalues # store for later

  bounded <- FALSE
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
  }else{
    opt.method <- misc.options$mono.method # FTP: New bit of info that must be
                                           # supplied through control
    if (!(opt.method %in% c("solnp", "slsqp"))) {
      stop("The optimiser method for contraint optimisation in R should be 'slsqp' or 'solnp'.")
    }
  }
  
  # if monotonicity has been requested but we are using key only then just
  # use optimx
  if(misc.options$mono & is.null(ddfobj$adjustment)){
    misc.options$mono <- FALSE
    # set back to default optimiser
    opt.method <- optim.options$optimx.method
  }

  # if nofit=TRUE, we just want to set the parameters, calculate the
  # likelihood and exit
  if(misc.options$nofit){
    if(!is.null(initialvalues) && any(is.na(initialvalues))){
      stop("No fitting, but initial values not specified!\n")
    }

    lt <- list()
    lt$par <- initialvalues
    lt$value <- flnl(initialvalues, ddfobj, misc.options)
    lt$hessian <- NULL
    lt$model <- list(scalemodel=misc.options$scalemodel)
    lt$converge <- 0
    lt$message <- "MAYBE CONVERGENCE?"

    lt$aux <- c(optim.options, bounds, misc.options)
    ddfobj <- assign.par(ddfobj, lt$par)
    lt$aux$ddfobj <- ddfobj

    lt$optim.history <- rbind(misc.options$optim.history,
                              c(lt$conv, -lt$value, lt$par))

    lt$conv <- NULL
    return(lt)
  }

  # save last value of the lnl -- starting value
  lnl.last <- Inf

  # recover the optimisation history
  optim.history <- misc.options$optim.history

  # save the list of optimisation methods
  opt.method.save <- opt.method

  # Continue fitting until convergence occurs while parameters are within
  # their bounds

  itconverged <- FALSE
  while(!itconverged){
    # Call optimization routine to find constrained mles; upon
    # completion add the user specified models and return the list.

    ### Monotonically constrained optimisation
    if(misc.options$mono){

      # fail if int.range is a matrix
      if(is.matrix(misc.options$int.range) && nrow(misc.options$int.range)>1){
        stop("Monotonicity constraints not available with multiple integration ranges")
      }

      # lower and upper bounds of the inequality constraints
      lowerbounds.ic <- rep(0, 2*misc.options$mono.points)
      upperbounds.ic <- rep(10^10, 2*misc.options$mono.points)

      # small initialvalues lead to errors in solnp, so work around that
      if (opt.method == "solnp") {
        initialvalues[initialvalues<1e-2] <- sign(initialvalues[initialvalues<1e-2]) * 1e-2
      }
      if (showit == 0) {
        ## The derivative-based SLSQP solver
        if (opt.method == "slsqp") {
          lt <- suppressWarnings(
            try(nloptr(x0 = initialvalues,
                       eval_f = flnl,
                       eval_g_ineq = flnl.constr.neg,
                       eval_grad_f = flnl.grad,
                       eval_jac_g_ineq = flnl.constr.grad.neg,
                       lb = lowerbounds, ub = upperbounds,
                       opts = list(ftol_rel = misc.options$mono.tol, 
                                   ftol_abs = 0.0, 
                                   xtol_rel = 0.0,
                                   maxeval = 1000,
                                   print_level = as.integer(showit),
                                   algorithm = "NLOPT_LD_SLSQP"),
                       ddfobj = ddfobj,
                       misc.options = misc.options,
                       fitting = "all"),
                silent = TRUE))
          if (!inherits(lt, "try-error")) {
            if (lt$status >= 0) { # https://jhelvy.github.io/logitr/reference/statusCodes.html
              lt$convergence <- 0 
            } else {
              lt$convergence <- 1 # Have 1 now as code for failed convergence, not sure if correct
              # opt.method <- "auglag"
              # cat("The SLSQP solver did not converge from starting values:", initialvalues, 
              #     "!\nTrying an augmented Lagrangian with interior SLSQP solver instead.\n")
            } 
          }
        }
        ## The derivative-free augmented Lagrangian solver by Yinyu Ye.
        if (opt.method == "solnp") {
          lt <- suppressWarnings(
            try(solnp(pars=initialvalues, fun=flnl, eqfun=NULL, eqB=NULL,
                      ineqfun=flnl.constr,
                      ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                      LB=lowerbounds, UB=upperbounds,
                      ddfobj=ddfobj, misc.options=misc.options,
                      control=list(trace=as.integer(showit),
                                   tol=misc.options$mono.tol,
                                   delta=misc.options$mono.delta,
                                   outer.iter=misc.options$mono.outer.iter)
            ),
            silent=TRUE))
        } 
        if (!(opt.method %in% c("solnp", "auglag", "slsqp"))) {
          stop("Constraint solver is not 'auglag', 'solnp', or 'slsqp'!")
        }
      } else {
        ## The derivative-based SLSQP solver
        if (opt.method == "slsqp") {
          lt <- try(nloptr(x0 = initialvalues, 
                       eval_f = flnl,
                       eval_g_ineq = flnl.constr.neg,
                       eval_grad_f = flnl.grad,
                       eval_jac_g_ineq = flnl.constr.grad.neg,
                       lb = lowerbounds, ub = upperbounds,
                       opts = list(ftol_rel = misc.options$mono.tol, 
                                   ftol_abs = 0.0, 
                                   xtol_rel = 0.0,
                                   maxeval = 1000,
                                   print_level = as.integer(showit),
                                   algorithm = "NLOPT_LD_SLSQP"),
                       ddfobj = ddfobj,
                       misc.options = misc.options,
                       fitting = "all"))
          if (!inherits(lt, "try-error")) {
            if (lt$status >= 0) { # https://jhelvy.github.io/logitr/reference/statusCodes.html
              lt$convergence <- 0 
            } else {
              lt$convergence <- 1 # Have 1 now as code for failed convergence, not sure if correct
              # opt.method <- "auglag"
              # warning("The SLSQP solver did not converge from starting values:", initialvalues, 
              #     "!\nTrying an augmented Lagrangian with interior SLSQP solver instead.\n")
            }
          }
        }
        ## The derivative-free augmented Lagrangian solver by Yinyu Ye
        if (opt.method == "solnp") {
          lt <- try(solnp(pars=initialvalues, fun=flnl, eqfun=NULL, eqB=NULL,
                          ineqfun=flnl.constr,
                          ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                          LB=lowerbounds, UB=upperbounds,
                          ddfobj=ddfobj, misc.options=misc.options,
                          control=list(trace=as.integer(showit),
                                       tol=misc.options$mono.tol,
                                       delta=misc.options$mono.delta,
                                       outer.iter=misc.options$mono.outer.iter)
          ))
        } 
        if (!(opt.method %in% c("solnp", "auglag", "slsqp"))) {
          stop("Constraint solver is not 'auglag', 'solnp', or 'slsqp'!")
        }
      }

      # only do something more complicated if we didn't converge above!
      if(inherits(lt, "try-error") || lt$convergence!=0 ){
        # we can use the gosolnp() function to explore the parameter space
        # randomly... Only use this if solnp was specified as constr solver
        if(misc.options$mono.random.start & opt.method == "solnp"){
          if(length(initialvalues)>1){
            # gosolnp doesn't work when there is only 1 parameter
            # since there is a bug that leaves the optimisation with a vector
            # when it is expecting a matrix. To avoid this bug, don't do the
            # multiple start points in that case (should only be unif+cos(1))

            if(showit==0){
              lt2 <- suppressWarnings(
                     try(gosolnp(pars=initialvalues, fun=flnl,
                                 eqfun=NULL, eqB=NULL,
                                 ineqfun=flnl.constr,
                                 ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                                 LB=lowerbounds, UB=upperbounds,
                                 ddfobj=ddfobj, misc.options=misc.options,
                                 control=list(trace=as.integer(showit),
                                              tol=misc.options$mono.tol,
                                              delta=misc.options$mono.delta,
                                              outer.iter=misc.options$mono.outer.iter),
                                 distr = rep(1, length(lowerbounds)),
                                 n.restarts = 2, n.sim = 200,
                                 rseed=as.integer(runif(1)*1e9)),
                          silent=TRUE))
            }else{
              lt2 <- try(gosolnp(pars=initialvalues, fun=flnl,
                                 eqfun=NULL, eqB=NULL,
                                 ineqfun=flnl.constr,
                                 ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                                 LB=lowerbounds, UB=upperbounds,
                                 ddfobj=ddfobj, misc.options=misc.options,
                                 control=list(trace=as.integer(showit),
                                              tol=misc.options$mono.tol,
                                              delta=misc.options$mono.delta,
                                              outer.iter=misc.options$mono.outer.iter),
                                 distr = rep(1, length(lowerbounds)),
                                 n.restarts = 2, n.sim = 200,
                                 rseed=as.integer(runif(1)*1e9)))
            }

            # was this better than the first time?
            if(!inherits(lt2, "try-error")){
              if(inherits(lt, "try-error") ||
                 (!is.na(lt2$values[length(lt2$values)]) &&
                 (lt2$values[length(lt2$values)] <
                  lt$values[length(lt$values)]))){
                lt <- lt2
              }
            } # end "was it better" check
          } # end par length check
        } else {
          # otherwise we do non-random par space exploration using a grid

          # n.sim=200 for gosolnp, so we want a comparable systematic grid
          # so we want n.notsim=200= (seq len)^(par dimensions)
          # => ~~ 200^(1/(par dimensions))
          # but still need some resolution on the grid, so make the min
          # sequence be length 4??
          n.notsim <- max(c(4, ceiling(200^(1/length(lowerbounds)))))
          # create a sequence in each parameters direction
          par_grid <- apply(cbind(lowerbounds, upperbounds), 1,
                            \(x) seq(x[1], x[2], length.out=n.notsim))
          # trim the extremes
          par_grid <- par_grid[-c(1, nrow(par_grid)), , drop=FALSE]
          # then expand out the options
          par_grid <- as.matrix(expand.grid(as.data.frame(par_grid)))
          # add in the initialvalues pars too
          par_grid <- rbind(par_grid, initialvalues)

          # small initialvalues lead to errors in solnp, so work around that
          if (opt.method == "solnp") {
            par_grid[abs(par_grid)<1e-2] <- sign(par_grid[abs(par_grid)<1e-2])*1e-2
          }

          flnl_wrap <- function(...){
            e <- try(flnl(...), silent=TRUE)
            if(inherits(e, "try-error")){
              return(NA)
            }
            e
          }
          
          grid_lnls <- apply(par_grid, 1, flnl_wrap,
                             ddfobj=ddfobj, misc.options=misc.options)
          
          ## Check which grid points meet the constraint. 
          ## Problem is that sometimes none of the grid points meet the constraint. What then?
          ## Commented out below for now, but it might be needed?
          constr.met <- function(...) {
            return(all(flnl.constr(...) > 0))
          }
          grid_constr_met <- apply(par_grid, 1, constr.met,
                                    ddfobj=ddfobj, misc.options=misc.options)
          
          ## If none of the grid points meet the constraint, save the most 
          ## recent values and break out of the fitting process. 
          if (sum(grid_constr_met) == 0) {
            # if we ran out of refits we need to get out of here
            # lt <- list()
            lt$conv <- 1
            lt$value <- lnl.last
            lt$par <- initialvalues
            optim.history <- rbind(optim.history, c(lt$conv, -lt$value, lt$par)) 
            lt$message <- "FALSE CONVERGENCE & UNABLE TO FIND ALTERNATIVE STARTING VALUES"
            
            if(showit >= 2){
              cat("DEBUG: Optimisation failed, breaking and returning the object for analysis.\n")
            }
            break
          }
          grid_lnls <- grid_lnls[grid_constr_met]
          
          ## Find the grid point with the lowest negative lnl
          igrid <- which.min(grid_lnls)
          
          ## Now run it at best values using the right solver.
          if (showit == 0) {
            ## The derivative-based SLSQP solver
            if (opt.method == "slsqp") {
              lt <- suppressWarnings(
                try(nloptr(x0 = par_grid[igrid, ], 
                           eval_f = flnl,
                           eval_g_ineq = flnl.constr.neg,
                           eval_grad_f = flnl.grad,
                           eval_jac_g_ineq = flnl.constr.grad.neg,
                           lb = lowerbounds, ub = upperbounds,
                           opts = list(ftol_rel = misc.options$mono.tol, 
                                       ftol_abs = 0.0, 
                                       xtol_rel = 0.0,
                                       maxeval = 1000,
                                       print_level = as.integer(showit),
                                       algorithm = "NLOPT_LD_SLSQP"),
                           ddfobj = ddfobj,
                           misc.options = misc.options,
                           fitting = "all"), 
                    silent = TRUE))
              if (!inherits(lt, "try-error")) {
                if (lt$status >= 0) { # https://jhelvy.github.io/logitr/reference/statusCodes.html
                  lt$convergence <- 0 
                } else {
                  lt$convergence <- 1 # Have 1 now as code for failed convergence, not sure if correct
                }
              }
            }
            if (opt.method == "solnp") {
              lt <- suppressWarnings(
                try(solnp(pars = par_grid[igrid, ], fun=flnl, eqfun=NULL, eqB=NULL,
                          ineqfun=flnl.constr,
                          ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                          LB=lowerbounds, UB=upperbounds,
                          ddfobj=ddfobj, misc.options=misc.options,
                          control=list(trace=as.integer(showit),
                                       tol=misc.options$mono.tol,
                                       delta=misc.options$mono.delta,
                                       outer.iter=misc.options$mono.outer.iter)
                ),
                silent=TRUE))
            } 
          } else {
            ## The derivative-based SLSQP solver
            if (opt.method == "slsqp") {
              lt <- try(nloptr(x0 = par_grid[igrid, ], 
                           eval_f = flnl,
                           eval_g_ineq = flnl.constr.neg,
                           eval_grad_f = flnl.grad,
                           eval_jac_g_ineq = flnl.constr.grad.neg,
                           lb = lowerbounds, ub = upperbounds,
                           opts = list(ftol_rel = misc.options$mono.tol, 
                                       ftol_abs = 0.0, 
                                       xtol_rel = 0.0,
                                       maxeval = 1000,
                                       print_level = as.integer(showit),
                                       algorithm = "NLOPT_LD_SLSQP"),
                           ddfobj = ddfobj,
                           misc.options = misc.options,
                           fitting = "all"))
              if (!inherits(lt, "try-error")) {
                if (lt$status >= 0) { # https://jhelvy.github.io/logitr/reference/statusCodes.html
                  lt$convergence <- 0 
                } else {
                  lt$convergence <- 1 # Have 1 now as code for failed convergence, not sure if correct
                }
              }
            }
            if (opt.method == "solnp") {
              lt <- try(solnp(pars=par_grid[igrid, ], fun=flnl, eqfun=NULL, eqB=NULL,
                              ineqfun=flnl.constr,
                              ineqLB=lowerbounds.ic, ineqUB=upperbounds.ic,
                              LB=lowerbounds, UB=upperbounds,
                              ddfobj=ddfobj, misc.options=misc.options,
                              control=list(trace=as.integer(showit),
                                           tol=misc.options$mono.tol,
                                           delta=misc.options$mono.delta,
                                           outer.iter=misc.options$mono.outer.iter)
              ))
            } 
          }
        } # end random vs. non-random par space exploration
      } # end if constraint solver didn't converge the first time

      # If that failed then make a dummy object
      if(inherits(lt, "try-error")){
        lt <- list()
        lt$conv <- 9
        lt$value <- lnl.last
        lt$par <- initialvalues

        if(showit >= 2){
          cat("DEBUG: Optimisation failed, ignoring and carrying on...\n")
        }
      }else{
        lt$conv <- lt$convergence
        lt$message <- ""
        
        if (opt.method == "solnp") {
          lt$par <- lt$pars
          lt$value <- lt$values[length(lt$values)]
        } else if (opt.method %in% c("slsqp", "auglag")) {
          lt$par <- lt$solution
          lt$value <- lt$objective
        } else stop("Invalid constrained solver")
        
      }
    ## end monotonically constrained estimation
    } else {
    ## unconstrained optimisation

      # get the list of optimisation methods
      opt.method <- opt.method.save

      # SANN is not an optimx method, so do that first using optim
      if(any(opt.method == "SANN")){
        lt <- try(optim(initialvalues, flnl, method="SANN",
                        control=optim.options,
                        hessian=TRUE, lower=lowerbounds,
                        upper=upperbounds, ddfobj=ddfobj, fitting=fitting,
                        misc.options=misc.options), silent=TRUE)
        opt.method <- opt.method[opt.method != "SANN"]
        if(!inherits(lt, "try-error")){
          initialvalues <- lt$par
        }
      }

      # if one of the methods was nlminb and we need to rescale (see above)
      # then we need to call nlminb separately to ensure that the scaling is
      # passed (optimx will phase out scaling according to JC Nash, email
      # to DLM, 4 Dec 2014)
      # this code taken from Distance2

      # only do this if we have covariates
      # AND we are fitting all of the parameters
      # AND parscale is TRUE
      # AND using nlminb
      # rescaling only happens if the scaling factor is above some threshold
      # see ?rescale_pars
      if(!is.null(ddfobj$scale$formula) && (ddfobj$scale$formula != ~1) &&
         !is.null(optim.options$parscale) &&
         all(optim.options$parscale) &&
          (opt.method=="nlminb") && (fitting=="all")){

        # get the rescaling
        if(is.logical(optim.options$parscale)){
          optim.options$parscale <- rescale_pars(initialvalues, ddfobj)
        }

        # run the optimiser
        lt <- try(nlminb_wrapper(par=initialvalues, ll=flnl,
                                 lower=lowerbounds,
                                 upper=upperbounds,
                                 mcontrol=optim.options, ddfobj=ddfobj,
                                 misc.options=misc.options))
        # remove nlminb from the list
        opt.method <- opt.method[opt.method != "nlminb"]

        # get new initialvalues
        if(inherits(lt, "try-error") || any(is.na(lt[,1:attr(lt,"npar")]))){
          if(showit >= 2){
            cat("DEBUG: Optimisation failed, ignoring and carrying on...\n")
          }
        }
        # put the results in a nice format
        lt <- parse.optimx(lt, lnl.last, initialvalues)
        initialvalues <- lt$par
        # reset parscale to be logical, it will be recalculated each time
        optim.options$parscale <- TRUE
      } # end rescaled nlminb

      ## use optimx if there are methods left
      if(length(opt.method)>0){

        # remove scaling, as otherwise there's weird conflicts
        optim.options$parscale <- NULL
        # remove this again, because of conflicts
        optim.options$optimx.method <- NULL

        # if we're only fitting one part of the model, we need to do the
        # optimisation over that bit only (with only the relevant parameters
        # and corresponding bounds)
        if(fitting == "key"){
          savedvalues <- initialvalues
          parind <- getpar(ddfobj, index=TRUE)
          initialvalues <- initialvalues[1:parind[2]]
          ub <- upperbounds[1:parind[2]]
          lb <- lowerbounds[1:parind[2]]
        }else if(fitting=="adjust"){
          savedvalues <- initialvalues
          parind <- getpar(ddfobj, index=TRUE)
          initialvalues <- initialvalues[(parind[2]+1):parind[3]]
          ub <- upperbounds[(parind[2]+1):parind[3]]
          lb <- lowerbounds[(parind[2]+1):parind[3]]
        }else{
          savedvalues <- initialvalues
          ub <- upperbounds
          lb <- lowerbounds
        }

        # now run the optimizer
        lt <- try(optimx(initialvalues, flnl, method=opt.method,
                         control=optim.options, hessian=TRUE,
                         lower=lb, upper=ub,
                         ddfobj=ddfobj, fitting=fitting,
                         misc.options=misc.options), silent=TRUE)

        if(inherits(lt, "try-error") || any(is.na(lt[, 1:attr(lt,"npar")]))){
          if(showit >= 2){
            cat("DEBUG: Optimisation failed, ignoring and carrying on...\n")
          }
        }

        # put the results in a nice format
        lt <- parse.optimx(lt, lnl.last, initialvalues)

        # ensure that ddfobj has the full parameter set by putting the
        # saved values back in place
        if(fitting == "key"){
          if(parind[1] == 0){
            # half-normal
            ddfobj$scale$parameters <- lt$par
          }else{
            # hazard-rate
            ddfobj$shape$parameters <- lt$par[1]
            ddfobj$scale$parameters <- lt$par[2:parind[2]]
          }
          if(parind[3] != 0){
            lt$par <- c(lt$par, savedvalues[(parind[2]+1):parind[3]])
          }
        }else if(fitting == "adjust"){
          ddfobj$adjustment$parameters <- lt$par
          lt$par <- c(savedvalues[1:parind[2]], lt$par)
        }

        initialvalues <- lt$par
      } # end optimx
    } # end unconstrained optimisation

    # now we're done with optimization, what happened?
    # Print debug information
    if(showit>=2){
      cat("DEBUG: Converge   =", lt$conv, "(", fitting, ")\n",
          "      nll        =", lt$value,"\n",
          "      parameters =", paste(round(lt$par, 7), collapse=", "), "\n")
    }

    if(is.null(lt$value)){
      lt$value <- NA
      lt$conv <- 2
    }
    optim.history <- rbind(optim.history, c(lt$conv, -lt$value, lt$par)) 

    # check whether parameters hit their bounds
    bounded <- check.bounds(lt, lowerbounds, upperbounds, ddfobj,
                            showit, setlower, setupper)

    # did the optimisation converge? (code 0 is GOOD)
    if(!bounded & (lt$conv==0 | !refit)){
      itconverged <- TRUE
    }else{
    # If we don't have convergence what do we do

      # first check to see if any parameters were NA and if we're in debug mode
      # return that bad object for diagnosis
      if(any(is.na(lt$par)) & misc.options$debug){
        # if there was no convergence then just return the
        # lt object for debugging
        warning("Problems with fitting model. Did not converge")
        lt$optim.history <- optim.history
        break
      }


      # increment the refit counter and check we've not run out of refits
      refit.count <- refit.count + 1
      if(refit.count < nrefits){
        if(showit >= 1){
          cat("DEBUG: No convergence. Refitting ...\n")
        }

        # move starting values around if they were not close to the boundary
        # (if they were at the upper/lower bounds then we should make the
        # bounds wider rather than mess with the initial values, see below)
        if(!bounded & fitting=="all"){
          initialvalues <- lt$par*runif(length(initialvalues), 
                                        sign(lowerbounds-1),
                                        sign(upperbounds+1))
          #initialvalues <- runif(length(initialvalues),
          #                              lowerbounds,
          #                              upperbounds)
        }
        lnl.last <- lt$value

        # if had NAs in lt$par, replace them with the values we started with
        initialvalues[is.na(initialvalues)] <-
                                    initialvalues.set[is.na(initialvalues)]
      }else{
        # if we ran out of refits we need to get out of here
        lt$message <- "FALSE CONVERGENCE"
        lt$conv <- 1
        break
      }
    } # end convergence check

    # for parameters that were at their bounds, we widen those bounds
    if(bounded){
      bound.low <- abs(lt$par-lowerbounds)<1e-6
      bound.hi <- abs(lt$par-upperbounds)<1e-6
      bound.scale <- 0.5
      if(!setlower) {
        lowerbounds[bound.low] <- lowerbounds[bound.low] -
                                   bound.scale*abs(lowerbounds[bound.low])
        lowerbounds[bound.low & lowerbounds>0 & lowerbounds < 0.5] <- -0.5
      }
      if(!setupper){
        upperbounds[bound.hi] <- upperbounds[bound.hi] +
                                  bound.scale*abs(upperbounds[bound.hi])
        upperbounds[bound.hi & upperbounds<0 & upperbounds > -0.5] <- 0.5
      }

      if(showit>=1){
        cat("DEBUG: Refitting ...\n")
      }
    } # end of bounds checking
  } # end of while(!bounded & !itconverged)

  # build the lt object to return
  lt$aux <- c(optim.options, bounds, misc.options)
  ddfobj <- assign.par(ddfobj, lt$par)
  lt$aux$ddfobj <- ddfobj
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
